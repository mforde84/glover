# glover - genelist_overlaps

Given a geneset list formatted according to mSigDB standards, output overlapping genesets for any number of specified pathways against the entire geneset list.

compile: g++ -o glover genelist_overlaps.cpp --std=c++17

usage: ./glover <mSigDB.gmt file OR similarily formated> <search string 1> <search string 2> ... > output.csv

eg, ./glover gene.lists ZHANG_INTERFERON_RESPONSE AACTGAC,MIR-223 > example.csv

