#!/usr/bin/env julia

using Random
using Distributions
using StatsBase
using DelimitedFiles
using Dates
using Parameters


using SQLite: DB, Stmt, bind!
using SQLite.DBInterface: execute
using JSON
