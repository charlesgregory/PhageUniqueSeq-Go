package main

import "os"
import (
	"PhageAnalysis"
)


func main() {
	PhageAnalysis.ParseArgs(os.Args)
	//PhageAnalysis.DoPrimerAnalysis(18,25)
	//PhageAnalysis.PrimerAnalysisTest(18,"Mycobacterium")
	//x:=PhageAnalysis.ReadUniquePrimers("A1","Mycobacterium")
	//println(len(x))
	//PhageAnalysis.MatchPrimers()
	PhageAnalysis.MatchingTest("Mycobacterium","A1")

}
