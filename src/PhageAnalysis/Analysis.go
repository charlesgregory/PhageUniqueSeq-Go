package PhageAnalysis

import (
	"bufio"
	//"time"
	//"fmt"
	//"strconv"
	//"math"
	"os"
	"log"
	"time"
	"strconv"
	"fmt"
)
/***/
func PrimerAnalysis(bps int, w *bufio.Writer,phageList map[string]map[string]map[string]string,para bool){
	//strain:="Mycobacterium"
	//clusters:=phageList[strain]
	var strain,cluster string
	var clusters map[string]map[string]string
	var phages map[string]string
	var primers map[uint64]Primer
	var temp map[uint64]Primer
	var clustersMap map[uint8]string
	for strain,clusters= range phageList{
		t1:=time.Now()
		primers= make(map[uint64]Primer)
		clustersMap=make(map[uint8]string)
		fmt.Println(strain)
		var clusterNum uint8=0
		var primer uint64
		var x,prim Primer
		var check,check2 bool
		var primChan =make(chan map[uint64]Primer,len(clusters))
		if(para){
			for cluster,phages=range clusters {
				clustersMap[clusterNum]=cluster
				go byCluster(phages,clusterNum,bps,clustersMap,cluster,primChan)
				clusterNum=clusterNum+1

			}
		}else{
			for cluster,phages=range clusters {
				clustersMap[clusterNum]=cluster
				byCluster(phages,clusterNum,bps,clustersMap,cluster,primChan)
				clusterNum=clusterNum+1

			}
		}
		for range clusters{
			temp=<-primChan
			for primer,prim=range temp{
				_, check = primers[primer]
				_, check2 = primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]
				if (!check&&!check2) {
					primers[primer] = prim
				} else {
					if (check) {
						x = primers[primer]
						x.combine(prim)
						primers[primer]=x
					}
					if(check2){
						x = primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]
						x.combine(prim)
						primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]=x
					}

				}
			}
		}
		t2:=time.Now()
		fmt.Println(time.Since(t1).Minutes())
		fmt.Println(len(primers))
		//keepPrimers:=make(map[uint64]Primer)
		count:=0
		for p,v:=range primers{
			if(len(v.clusters)==1){
				//count++
				var primerClust uint8
				for key,_:=range v.clusters{
					primerClust=key
				}
				if(len(clusters[clustersMap[primerClust]])==int(v.phagecount)){
					count++
					w.WriteString(strain+",")
					w.WriteString(clustersMap[primerClust])
					w.WriteString(",")
					w.WriteString(strconv.FormatUint(p,10)+"\n")
				}
			}
		}
		fmt.Println(time.Since(t2).Minutes())
		fmt.Println(count)


	}

}
func byCluster(phages map[string]string,
clusterNum uint8,bps int,clustersMap map[uint8]string,
cluster string, primerChan chan map[uint64]Primer){
	fmt.Println("\t"+cluster)
	var seq string
	var check,check2 bool
	var phagprimers map[uint64]bool
	var primers =make(map[uint64]Primer)
	var primer uint64
	//var prev uint16=0
	for _,seq=range phages{
		//fmt.Println("\t\t"+phage)
		//go makePrimers(seq,bps,primChan)
		phagprimers=make(map[uint64]bool)
		for k:= 0; k <= len(seq) - bps; k++ {
			primer=twoBitEncoding(seq[k:k+bps])
			_, check = phagprimers[primer]
			_, check2 = phagprimers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]
			if (!check&&!check2) {
				phagprimers[primer]=true
			}
		}
		for primer, _ = range phagprimers {
			_, check = primers[primer]
			_, check2 = primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]
			if (!check&&!check2) {
				var x Primer
				x.addCluster(clusterNum)
				x.addPhage()
				primers[primer] = x
			} else {
				if (check) {
					var x Primer;
					x = primers[primer]
					x.addPhage()
					primers[primer]=x
				}else if(check2){
					var x Primer;
					x = primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]
					x.addPhage()
					primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]=x
				}

			}
			//if(twoBitDecode(primer)=="CTTCCACGGCGAGGACCC"||
			//	RevComplement(twoBitDecode(primer))=="CTTCCACGGCGAGGACCC"){
			//	if(primers[primer].phagecount-prev>1){
			//		println()
			//	}
			//	prev=primers[primer].phagecount
			//}
			//if(cluster=="A1"&&(twoBitDecode(primer)=="TGAGAGCCCCGTAGACGG")){
			//	println(primers[primer].phagecount)
			//}

		}
		//fmt.Println(phage)
	}
	primerChan<-primers
}
func PrimerAnalysisWorker(bps int, w *bufio.Writer,phageList map[string]map[string]map[string]string,threads int){
	//strain:="Mycobacterium"
	//clusters:=phageList[strain]
	var strain,cluster string
	var clusters map[string]map[string]string
	var phages map[string]string
	var primers map[uint64]Primer
	var temp map[uint64]Primer
	var clustersMap map[uint8]string
	for strain,clusters= range phageList{
		t1:=time.Now()
		primers= make(map[uint64]Primer)
		clustersMap=make(map[uint8]string)
		fmt.Println(strain)
		var clusterNum uint8=0
		var primer uint64
		var x,prim Primer
		var check,check2 bool
		var jobs = make(chan AnalysisWork, len(clusters))
		var primChan =make(chan map[uint64]Primer,len(clusters))


		for w := 1; w <= threads; w++ {
			go AnalysisWorker(jobs, primChan)
		}

		for cluster,phages=range clusters {
			clustersMap[clusterNum]=cluster
			jobs <- AnalysisWork{phages,clusterNum,bps,cluster}
			clusterNum=clusterNum+1

		}
		close(jobs)
		for range clusters{
			temp=<-primChan
			for primer,prim=range temp{
				_, check = primers[primer]
				_, check2 = primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]
				if (!check&&!check2) {
					primers[primer] = prim
				} else {
					if (check) {
						x = primers[primer]
						x.combine(prim)
						primers[primer]=x
					}
					if(check2){
						x = primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]
						x.combine(prim)
						primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]=x
					}

				}
			}
		}
		t2:=time.Now()
		fmt.Println(time.Since(t1).Minutes())
		fmt.Println(len(primers))
		//keepPrimers:=make(map[uint64]Primer)
		count:=0
		for p,v:=range primers{
			if(len(v.clusters)==1){
				//count++
				var primerClust uint8
				for key,_:=range v.clusters{
					primerClust=key
				}
				if(len(clusters[clustersMap[primerClust]])==int(v.phagecount)){
					count++
					w.WriteString(strain+",")
					w.WriteString(clustersMap[primerClust])
					w.WriteString(",")
					w.WriteString(strconv.FormatUint(p,10)+"\n")
				}
			}
		}
		fmt.Println(time.Since(t2).Minutes())
		fmt.Println(count)


	}

}
func DoPrimerAnalysis(from int, to int, para bool){
	f, err := os.Create(WorkingDir +"Data"+pathslash+"Unique.csv")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	phageList:=ParsePhages()
	for i:=from;i<=to;i++{
		PrimerAnalysis(i,w,phageList,para)
	}
	err = w.Flush() // Don't forget to flush!
	if err != nil {
		log.Fatal(err)
	}
}
func TestPrimerAnalysis(from int, to int,threads int){
	f, err := os.Create(WorkingDir +"Data"+pathslash+"Unique.csv")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	phageList:=ParsePhages()
	for i:=from;i<=to;i++{
		PrimerAnalysisWorker(i,w,phageList,threads-1)
	}
	err = w.Flush() // Don't forget to flush!
	if err != nil {
		log.Fatal(err)
	}
}