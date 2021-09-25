dbSimInfo = SQLite.DB(dbSimInfoPath)
tableNamesTypes = ["$(table_info.name) $(table_info.type)" for table_info in execute(dbSimInfo,"PRAGMA table_info(param_combos)")]
tableNamesTypes = join(tableNamesTypes,", ")

execute(dbOutput, "CREATE TABLE runs (run_id INTEGER, combo_id INTEGER, replicate INTEGER)")
execute(dbOutput, "CREATE TABLE param_combos ($(tableNamesTypes...))")

tableNames = ["$(table_info.name)" for table_info in execute(dbSimInfo,"PRAGMA table_info(param_combos)")]
tableNames = join(tableNames,", ")
execute(dbOutput, "BEGIN TRANSACTION")
execute(dbOutput,"ATTACH DATABASE '$(dbSimInfoPath)' as dbSimInfo")
execute(dbOutput,"INSERT INTO param_combos($(tableNames)) SELECT * FROM dbSimInfo.param_combos")
execute(dbOutput,"INSERT INTO runs (run_id, combo_id, replicate) SELECT run_id, combo_id, replicate FROM dbSimInfo.runs")
execute(dbOutput, "COMMIT")
