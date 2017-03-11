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
	//PhageAnalysis.DoPrimerAnalysis(18,25,true)
	//PhageAnalysis.PrimerAnalysisTest(18,"Mycobacterium")
	//x:=PhageAnalysis.ReadUniquePrimers("B1","Mycobacterium")
	//println(len(x))
	//PhageAnalysis.MatchPrimersParallel()
	//str:=PhageAnalysis.ReadFile(PhageAnalysis.WorkingDir+"\\Fastas\\20ES.fasta")
	//println(str)
	//PhageAnalysis.MatchingTest("Mycobacterium","A1")
	PhageAnalysis.UniqueConfirm("Mycobacterium")
}
