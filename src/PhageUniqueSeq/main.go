package main

import "os"
import (
	"PhageAnalysis"
	"bufio"
	//"strings"
	//"strconv"
	//"fmt"
	//"strings"
	//"strconv"
	"regexp"
	//"fmt"
)


func main() {
	//PhageAnalysis.DownloadFromUrl("http://phagesdb.org/api/sequenced_phages/",
	//	"C:\\Users\\musta_000\\IdeaProjects\\GoProjects\\Fastas\\Phagelist.txt")


	file,_:=os.Open("C:\\Users\\musta_000\\IdeaProjects\\GoProjects\\Fastas\\Phagelist.txt")
	scanner := bufio.NewReader(file)
	line,_,_:=scanner.ReadLine()
	str:=string(line)
	match, _ := regexp.Compile("\"count\":....")
	for _,x:=range match.FindAllString(str,-1){
		println(x)
	}
	match, _ = regexp.Compile("\"phage_name\":\".*\",\"o")
	for _,x:=range match.FindAllString(str,-1){
		println(x)
	}
	//count,_:=strconv.ParseInt(str[strings.Index(str,"\"count\"")+len("\"count\"")+1:strings.Index(str,",")],10,32)
	//next:=str[strings.Index(str,"\"next\"")+len("\"next\"")+2:strings.Index(str,",\"previous")-1]
	//println(count)
	//println(next)
	println(str)
	line,_,_=scanner.ReadLine()
	println(string(line))
	PhageAnalysis.ParseArgs(os.Args)
	//primers:=createClusterPrimerMap()
	//for strain,clusters:=range primers{
	//	println(strain)
	//	for cluster,arr:=range clusters{
	//		println("\t"+cluster)
	//		for i:=0;i<len(arr);i++{
	//			print("\t\t")
	//			println(arr[i])
	//		}
	//	}
	//}
	//fastas:=PhageAnalysis.GetFastas()
	//for strain,clusters:=range fastas{
	//	println(strain)
	//	for cluster,phages:=range clusters{
	//		println("\t"+cluster)
	//		for phage,seq:=range phages{
	//			println("\t\t"+phage)
	//			print("\t\t")
	//			println(len(seq))
	//		}
	//	}
	//}
	//doPrimerAnalysis()
	//matchPrimers()
	//primerAnalysisTest()
	//println(twoBitDecode(uint64(788971462002)))
	//for k,v:=range readUniqueClusters(){
	//	println(k)
	//	println(len(v))
	//	for j,_:=range v{
	//		println(j)
	//	}
	//}
	//openUniqueCsv("percival")
	//DoAll()

}
