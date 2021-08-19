#!/usr/bin/env julia

using Pkg

#Pkg.add("JSON2")
Pkg.add("JSON")
Pkg.add("Parameters")
Pkg.add("StatsBase")
Pkg.add("Distributions")
Pkg.add("SQLite")
Pkg.add("Parameters")

using Random
using Distributions
using StatsBase
using DelimitedFiles
using Dates
using Parameters


using SQLite: DB, Stmt, bind!
using SQLite.DBInterface: execute
