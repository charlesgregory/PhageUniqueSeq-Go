package main

import "os"
import (
	"PhageAnalysis"
)


func main() {
	PhageAnalysis.ParseArgs(os.Args)
	//PhageAnalysis.TestRegex()
	//PhageAnalysis.DoPrimerAnalysis(18,18)
	//PhageAnalysis.PrimerAnalysisTest(18,"Mycobacterium")
	//x:=PhageAnalysis.ReadUniquePrimers("A1","Mycobacterium")
	//println(len(x))
	//PhageAnalysis.MatchPrimers()
	//str:=PhageAnalysis.ReadFile(PhageAnalysis.WorkingDir+"\\Fastas\\20ES.fasta")
	//println(str)
	//PhageAnalysis.MatchingTest("Mycobacterium","A1")
	PhageAnalysis.UniqueConfirm("Mycobacterium")
}
