function fitness(allele,centerallele,maxallele,maxfitness,evofunction)
    if evofunction == 1
        minallele = 2*centerallele - maxallele
        scale = -1*maxfitness/((centerallele - minallele)*(centerallele - maxallele))
        return -1*(allele - minallele)*(allele - maxallele)*scale
    end
end
