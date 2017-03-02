package PhageAnalysis

import (
	"fmt"
	"os"
	//"log"
	"bufio"
	//"time"
	"strconv"
	"strings"
)

//func DoAll(){
//	//phageList:=createPhageMap()
//	for strain,clusters:= range phageList {
//		fmt.Println(strain)
//		for cluster, phages := range clusters {
//			fmt.Println("\t" + cluster)
//			c := make(chan map[uint64]bool)
//			primers :=make(map[uint64]bool)
//			for _, seq := range phages {
//				go makePrimers(seq,18,c)
//			}
//			for i := 0; i < len(phages); i++{
//				temp:=<-c
//				for k,v :=range temp{
//					primers[k]=v
//				}
//			}
//		}
//	}
//}

//func primerAnalysisTest(){
//	f, err := os.Create(path+"\\Data\\Test.csv")
//	if err != nil {
//		log.Fatal(err)
//	}
//	bps:=18
//	defer f.Close()
//	w := bufio.NewWriter(f)
//	phageList:=createPhageMap()
//	strain:="Mycobacterium"
//	//clusters:=phageList[strain]
//	cluster:="A1"
//	t1:=time.Now()
//	primers:= make(map[uint64]Primer)
//	clustersMap:=make(map[uint8]string)
//	fmt.Println(strain)
//	var clusterNum uint8=0
//	phages:=phageList[strain][cluster]
//	clustersMap[clusterNum]=cluster
//	fmt.Println("\t"+cluster)
//	//primChan := make(chan map[uint64]bool)
//	phagprimers:=make(map[uint64]bool)
//	for _,seq:=range phages{
//		//println("\t\t"+phage)
//		for k:= 0; k <= len(seq) - bps; k++ {
//			if(len(seq[k:k+bps])!=18){
//				println("cut")
//			}
//			primer:=twoBitEncoding(seq[k:k+bps])
//			if(len(twoBitDecode(primer))!=18&&len(seq[k:k+bps])==18){
//				println("decode")
//			}
//			phagprimers[primer]=true
//
//		}
//	}
//	//for i := 0; i < len(phages); i++ {
//	//	temp:=<-primChan
//	//	for k,v :=range temp{
//	//		phagprimers[k]=v
//	//	}
//	//}
//	for primer, _ := range phagprimers {
//		_, check := primers[primer]
//		_, check2 := primers[twoBitEncoding(revComplement(twoBitDecode(primer)))]
//		if (!check&&!check2) {
//			var x = make([]uint8, 1)
//			x[0] = clusterNum
//			primers[primer] = Primer{x, 1}
//		} else {
//			var x Primer;
//			if (check) {
//				x = primers[primer]
//			} else {
//				x = primers[twoBitEncoding(revComplement(twoBitDecode(primer)))]
//			}
//			x.phagecount=x.phagecount+1
//			var found bool = false
//			for i := 0; i < len(x.clusters); i++ {
//				if (clustersMap[x.clusters[i]] == cluster) {
//					found = true
//				}
//			}
//			if (!found) {
//				var newArr = make([]uint8, len(x.clusters) + 1)
//				for i := 0; i < len(x.clusters); i++ {
//					newArr[i] = x.clusters[i]
//				}
//				newArr[len(x.clusters)] = clusterNum
//				x.clusters=newArr
//			}
//			if (check) {
//				primers[primer]=x
//			} else {
//				primers[twoBitEncoding(revComplement(twoBitDecode(primer)))]=x
//			}
//
//		}
//
//	}
//	clusterNum=clusterNum+1
//	t2:=time.Now()
//	println(time.Since(t1).Minutes())
//	println(len(primers))
//	//keepPrimers:=make(map[uint64]Primer)
//	count:=0
//	for p,v:=range primers{
//
//		if(len(v.clusters)==1){
//			//count++
//			primerClust:=v.clusters[0]
//			//if(len(clusters[clustersMap[primerClust]])==v.phagecount){
//			count++
//			w.WriteString(strain+",")
//			w.WriteString(clustersMap[primerClust])
//			w.WriteString(",")
//			w.WriteString(strconv.FormatUint(p,10)+"\n")
//			//}
//		}
//	}
//	println(time.Since(t2).Minutes())
//	println(count)
//	err = w.Flush() // Don't forget to flush!
//	if err != nil {
//		log.Fatal(err)
//	}
//
//}
func readTest()[]uint64{
	f, _ := os.Open(workingDir +"\\Data\\Test.csv")
	var lines []uint64
	scanner := bufio.NewScanner(f)
	count:=0
	for scanner.Scan() {
		line:=scanner.Text()
		count++
		//if(strings.Contains(line,strain)&&strings.Contains(line,cluster)){
		strA:=strings.Split(line,",")
		i,_:=strconv.ParseUint(strA[2],10,64)
		lines = append(lines,i )
		//if(len(strA[2])!=18){
		//	println(strain+" "+cluster+" "+strA[2])
		//}
		//}

	}
	if err := scanner.Err(); err != nil {
		fmt.Fprintln(os.Stderr, err)
	}
	return lines
}