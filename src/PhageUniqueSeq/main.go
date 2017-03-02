package main

import "os"
import (
	"PhageAnalysis"
)


func main() {
	PhageAnalysis.ParseArgs(os.Args)
	PhageAnalysis.DoPrimerAnalysis()
	//x:=PhageAnalysis.ReadUniquePrimers("","Mycobacterium")
	//println(len(x))


}
