# This function counts number of peaks, logs time series of peaks their respective durations
function peakwallCount(uPercent, lPercent, shannonThreshold, dbOutput, dbTempSim)
    # Save as data frame
    series = DataFrame(execute(dbTempSim, "SELECT t, microbial_abundance FROM summary"))

    upperThreshold = floor(cCapacity * uPercent)
    lowerThreshold = floor(cCapacity * lPercent)

    thresholdVals = DataFrame(carrying_capacity=cCapacity,
        upper_percent=uPercent, lower_percent=lPercent,
        upper_threshold=upperThreshold, lower_threshold=lowerThreshold,
        shannon_threshold=shannonThreshold)
    thresholdVals |> SQLite.load!(dbOutput, "threshold_values", ifnotexists=true)

    # Identify crossings of thresholds
    coarseSeries = series[(series.microbial_abundance.>=upperThreshold).|(series.microbial_abundance.<=lowerThreshold), :]
    upper = 2 * (series.microbial_abundance .>= upperThreshold) # 2 signifies values above upper threshold
    lower = series.microbial_abundance .<= lowerThreshold # 1 signifies values below lower threshold, 0 for values between lower and upper threshold,
    coarseSeries = upper + lower

    # Initialize each table of output database as data frame
    peakSeries = DataFrame(t=series.t, peak_presence=Array{Int64,1}(zeros(length(series.t))))
    peakDurations = DataFrame(begin_row=Int64[], end_row=Int64[], begin_t=Float64[], end_t=Float64[],
        peak_number=Int64[], wall_number=Int64[], duration=Float64[], max_single_spacer_match_shan_div=Float64[])
    wallSeries = DataFrame(t=series.t, wall_presence=Array{Int64,1}(zeros(length(series.t))), bshannon=Array{Float64,1}(zeros(length(series.t))), vshannon=Array{Float64,1}(zeros(length(series.t))))

    peaks = 0
    walls = 0
    check = 0
    j1, j2, j3, j4 = 0, 0, 0, 0

    for i in 1:length(coarseSeries)
        println("Scanning through time series point $(i)")
        if i == length(coarseSeries)
            return peaks, walls, peakSeries, peakDurations, wallSeries
        end
        if coarseSeries[i] == 1 && coarseSeries[i+1] == 0
            check = 1
            j1 = i
        end
        if coarseSeries[i] == 0 && coarseSeries[i+1] == 2
            if check == 1
                check = 2
                j2 = i
            end
        end
        if coarseSeries[i] == 2 && coarseSeries[i+1] == 0
            if check == 2
                check = 3
                j3 = i
            end
        end
        if coarseSeries[i] == 0 && coarseSeries[i+1] == 1
            if check == 3
                j4 = i
                peaks = peaks + 1

                peakSeries.peak_presence[j1:j4] = Array{Int64,1}(ones(length(j1:j4)))
                peak_duration = peakSeries.t[j4] - peakSeries.t[j1]

                timeStmt = join(map(x -> string(" $(x)"), peakSeries.t[(j1):j4]),',')
                # timeStmt = append!([string("t = $(peakSeries.t[j1])")], timeStmt)
                bshannon = [shannon.bshannon for shannon in execute(dbTempShan, "SELECT bshannon FROM shannon_diversity WHERE t in ($(timeStmt...)) ORDER BY t")]
                vshannon = [shannon.vshannon for shannon in execute(dbTempShan, "SELECT vshannon FROM shannon_diversity WHERE t in ($(timeStmt...)) ORDER BY t")]
                spacerDiv = [div for (div,) in execute(
                    dbTempTri,
                    "SELECT shannon_diversity FROM single_match_diversity
                    WHERE t > $(peakSeries.t[j1]) AND t < $(peakSeries.t[j4])")]
                if maximum(bshannon) >= shannonThreshold
                    walls = walls + 1
                    wallSeries.wall_presence[j1:j4] = Array{Int64,1}(ones(length(j1:j4)))
                    wallSeries.bshannon[j1:j4] = bshannon[:, :]
                    wallSeries.vshannon[j1:j4] = vshannon[:, :]
                    if length(spacerDiv) > 0
                        push!(peakDurations, [j1 j4 peakSeries.t[j1] peakSeries.t[j4] peaks walls peak_duration maximum(spacerDiv)])
                    else
                        push!(peakDurations, [j1 j4 peakSeries.t[j1] peakSeries.t[j4] peaks walls peak_duration 0])
                    end
                else
                    if length(spacerDiv) > 0
                        push!(peakDurations, [j1 j4 peakSeries.t[j1] peakSeries.t[j4] peaks 0 peak_duration maximum(spacerDiv)])
                    else
                        push!(peakDurations, [j1 j4 peakSeries.t[j1] peakSeries.t[j4] peaks 0 peak_duration 0])
                    end
                end

                check = 0
            end
        end
    end
    return peaks, walls, peakSeries, peakDurations, wallSeries
    #end
end


function diversityToExtinction()
    (extV,) = execute(dbOutput, "SELECT viruses FROM extinction_occurrence")
    (walls,) = execute(dbOutput, "SELECT num_walls FROM microbial_peakwall_count")
    if extV[1] == 1 && walls[1] > 0
        wallNum = maximum([wallNum for (wallNum,) 
            in execute(dbOutput, 
                    "SELECT wall_number 
                    FROM microbial_peakwall_durations")])
        (timeCutOff,) = execute(dbOutput, 
                            "SELECT end_t 
                            FROM microbial_peakwall_durations
                            WHERE wall_number = $(wallNum)")
        timeCutOff = timeCutOff[1]
        maxDiv = maximum([div for (div,) 
                            in execute(dbTempTri, 
                            "SELECT shannon_diversity 
                            FROM single_spacer_match_diversity
                            WHERE t > $(timeCutOff)")])
        DataFrame(wall_viral_extinction=Int64[1], extinction_wall_number = Int64[wallNum+1], 
            max_single_spacer_match_shan_div=Float64[maxDiv]) |> SQLite.load!(dbOutput, "last_wall_extinction", ifnotexists=true)
    else
        DataFrame(wall_viral_extinction=Int64[0], extinction_wall_number = Int64[0], 
            max_single_spacer_match_shan_div=Float64[0]) |> SQLite.load!(dbOutput, "last_wall_extinction", ifnotexists=true)
    end
end