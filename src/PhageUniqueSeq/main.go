package main

import "os"
import (
	"PhageAnalysis"
	"time"
)


func main() {
	PhageAnalysis.ParseArgs(os.Args)
	//println(PhageAnalysis.WorkingDir)
	//PhageAnalysis.TestRegex()
	//PhageAnalysis.ExportClusterSummary()
	//PhageAnalysis.TestPrimer()
	t:=time.Now()
	PhageAnalysis.DoPrimerAnalysis(18,18,true)
	print("Time:")
	println(time.Since(t))
	//PhageAnalysis.PrimerAnalysisTest(18,"Mycobacterium")
	t=time.Now()
	PhageAnalysis.TestPrimerAnalysis(18,18,2)
	print("Time:")
	println(time.Since(t))
	//PhageAnalysis.MatchPrimersParallel()
	//str:=PhageAnalysis.ReadFile(PhageAnalysis.WorkingDir+"\\Fastas\\20ES.fasta")
	//println(str)
	//PhageAnalysis.MatchingTest("Mycobacterium","A1")
	//PhageAnalysis.UniqueConfirm("Mycobacterium")
	//PhageAnalysis.UniqueConfirmCluster("Mycobacterium","A1")
	//PhageAnalysis.FindPrimer("Mycobacterium","CTTCCACGGCGAGGACCC")
}
