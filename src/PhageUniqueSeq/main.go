package main

import "os"
import (
	"PhageAnalysis"
)


func main() {
	PhageAnalysis.ParseArgs(os.Args)
	//println(PhageAnalysis.WorkingDir)
	//PhageAnalysis.TestRegex()
	//PhageAnalysis.ExportClusterSummary()
	//PhageAnalysis.TestPrimer()
	//PhageAnalysis.DoPrimerAnalysis(18,25,true)
	//PhageAnalysis.PrimerAnalysisTest(18,"Mycobacterium")
	x:=PhageAnalysis.ReadMatchedPrimers("A1","Mycobacterium")
	var k PhageAnalysis.PrimerMatch
	var v PhageAnalysis.Stat
	for k,v=range x{
		print(k.F)
		print(",")
		print(k.R)
		print(",")
		print(v.Mean)
		print(",")
		print(v.Stddev)
		println()

	}
	//PhageAnalysis.MatchPrimersParallel()
	//str:=PhageAnalysis.ReadFile(PhageAnalysis.WorkingDir+"\\Fastas\\20ES.fasta")
	//println(str)
	//PhageAnalysis.MatchingTest("Mycobacterium","A1")
	//PhageAnalysis.UniqueConfirm("Mycobacterium")
	//PhageAnalysis.UniqueConfirmCluster("Mycobacterium","A1")
	//PhageAnalysis.FindPrimer("Mycobacterium","CTTCCACGGCGAGGACCC")
}
