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
type Primer struct{
	clusters []uint8
	phagecount int
}
type PrimerMatch struct {
	f uint64
	r uint64
}
func primerAnalysis(bps int, w *bufio.Writer){
	phageList:=ParsePhages()
	//strain:="Mycobacterium"
	//clusters:=phageList[strain]
	var strain,cluster,seq string
	var clusters map[string]map[string]string
	var phages map[string]string
	var primers map[uint64]Primer
	var phagprimers map[uint64]bool
	var clustersMap map[uint8]string
	//var primChan chan map[uint64]bool
	var primer uint64
	var check,check2 bool
	for strain,clusters= range phageList{
		t1:=time.Now()
		primers= make(map[uint64]Primer)
		clustersMap=make(map[uint8]string)
		fmt.Println(strain)
		var clusterNum uint8=0
		for cluster,phages=range clusters {
			clustersMap[clusterNum]=cluster
			fmt.Println("\t"+cluster)
			//primChan = make(chan map[uint64]bool)
			for _,seq=range phages{
				//fmt.Println("\t\t"+phage)
				//go makePrimers(seq,bps,primChan)
				phagprimers=make(map[uint64]bool)
				for k:= 0; k <= len(seq) - bps; k++ {
					primer=twoBitEncoding(seq[k:k+bps])
					phagprimers[primer]=true
				}
				for primer, _ = range phagprimers {
					_, check = primers[primer]
					_, check2 = primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]
					if (!check&&!check2) {
						var x = make([]uint8, 1)
						x[0] = clusterNum
						primers[primer] = Primer{x, 1}
					} else {
						var x Primer;
						if (check) {
							x = primers[primer]
						} else {
							x = primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]
						}
						x.phagecount=x.phagecount+1
						var found bool = false
						var clust uint8
						for _,clust=range x.clusters {
							if (clustersMap[clust] == cluster) {
								found = true
							}
						}
						if (!found) {
							var newArr = make([]uint8, len(x.clusters) + 1)
							for i := 0; i < len(x.clusters); i++ {
								newArr[i] = x.clusters[i]
							}
							newArr[len(x.clusters)] = clusterNum
							x.clusters=newArr
						}
						if (check) {
							primers[primer]=x
						} else {
							primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]=x
						}

					}

				}
				//fmt.Println(phage)
			}
			//for i := 0; i < len(phages); i++ {
			//	temp:=<-primChan
			//	fmt.Print(phage+" ")
			//	fmt.Println(len(temp))
			//	for k,v :=range temp{
			//		phagprimers[k]=v
			//	}
			//}
			//fmt.Println(len(phagprimers))
			clusterNum=clusterNum+1

		}
		t2:=time.Now()
		fmt.Println(time.Since(t1).Minutes())
		fmt.Println(len(primers))
		//keepPrimers:=make(map[uint64]Primer)
		count:=0
		for p,v:=range primers{

			if(len(v.clusters)==1){
				//count++
				primerClust:=v.clusters[0]
				if(len(clusters[clustersMap[primerClust]])==v.phagecount){
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
func primerAnalysisDouble(bps int, w *bufio.Writer){
	phageList:=ParsePhages()
	//strain:="Mycobacterium"
	//clusters:=phageList[strain]
	var strain,cluster,seq string
	var clusters map[string]map[string]string
	var phages map[string]string
	var primers map[uint64]Primer
	var phagprimers map[uint64]bool
	var clustersMap map[uint8]string
	//var primChan chan map[uint64]bool
	var primer,rprimer uint64
	var check bool
	for strain,clusters= range phageList{
		t1:=time.Now()
		primers= make(map[uint64]Primer)
		clustersMap=make(map[uint8]string)
		fmt.Println(strain)
		var clusterNum uint8=0
		for cluster,phages=range clusters {
			clustersMap[clusterNum]=cluster
			fmt.Println("\t"+cluster)
			//primChan = make(chan map[uint64]bool)
			for _,seq=range phages{
				//fmt.Println("\t\t"+phage)
				//go makePrimers(seq,bps,primChan)
				phagprimers=make(map[uint64]bool)
				for k:= 0; k <= len(seq) - bps; k++ {
					primer=twoBitEncoding(seq[k:k+bps])
					rprimer=twoBitEncoding(RevComplement(seq[k:k+bps]))
					phagprimers[primer]=true
					phagprimers[rprimer]=true
				}
				for primer, _ = range phagprimers {
					_, check = primers[primer]
					if (!check) {
						var x = make([]uint8, 1)
						x[0] = clusterNum
						primers[primer] = Primer{x, 1}
					} else {
						var x Primer;
						x = primers[primer]
						x.phagecount=x.phagecount+1
						var found bool = false
						var clust uint8
						for _,clust=range x.clusters {
							if (clustersMap[clust] == cluster) {
								found = true
							}
						}
						if (!found) {
							var newArr = make([]uint8, len(x.clusters) + 1)
							for i := 0; i < len(x.clusters); i++ {
								newArr[i] = x.clusters[i]
							}
							newArr[len(x.clusters)] = clusterNum
							x.clusters=newArr
						}
						primers[primer]=x

					}

				}
				//fmt.Println(phage)
			}
			//for i := 0; i < len(phages); i++ {
			//	temp:=<-primChan
			//	fmt.Print(phage+" ")
			//	fmt.Println(len(temp))
			//	for k,v :=range temp{
			//		phagprimers[k]=v
			//	}
			//}
			//fmt.Println(len(phagprimers))
			clusterNum=clusterNum+1

		}
		t2:=time.Now()
		fmt.Println(time.Since(t1).Minutes())
		fmt.Println(len(primers))
		//keepPrimers:=make(map[uint64]Primer)
		count:=0
		for p,v:=range primers{

			if(len(v.clusters)==1){
				//count++
				primerClust:=v.clusters[0]
				if(len(clusters[clustersMap[primerClust]])==v.phagecount){
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
func PrimerAnalysisParallel(bps int, w *bufio.Writer){
	phageList:=ParsePhages()
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
		for cluster,phages=range clusters {
			clustersMap[clusterNum]=cluster
			go byCluster(phages,clusterNum,bps,clustersMap,cluster,primChan)
			clusterNum=clusterNum+1

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
						x.phagecount=x.phagecount+prim.phagecount
					}
					if(check2){
						x = primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]
						x.phagecount=x.phagecount+prim.phagecount
					}
					var found bool = false
					var clust uint8
					for _,clust=range x.clusters {
						if (clust == prim.clusters[0]) {
							found = true
						}
					}
					if (!found) {
						var newArr = make([]uint8, len(x.clusters) + 1)
						for i := 0; i < len(x.clusters); i++ {
							newArr[i] = x.clusters[i]
						}
						newArr[len(x.clusters)] = prim.clusters[0]
						x.clusters=newArr
					}
					if (check) {
						primers[primer]=x
					}
					if(check2){
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
				primerClust:=v.clusters[0]
				if(len(clusters[clustersMap[primerClust]])==v.phagecount){
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
	for _,seq=range phages{
		//fmt.Println("\t\t"+phage)
		//go makePrimers(seq,bps,primChan)
		phagprimers=make(map[uint64]bool)
		for k:= 0; k <= len(seq) - bps; k++ {
			primer=twoBitEncoding(seq[k:k+bps])
			phagprimers[primer]=true
		}
		for primer, _ = range phagprimers {
			_, check = primers[primer]
			_, check2 = primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]
			if (!check&&!check2) {
				var x = make([]uint8, 1)
				x[0] = clusterNum
				primers[primer] = Primer{x, 1}
			} else {
				var x Primer;
				if (check) {
					x = primers[primer]
					x.phagecount=x.phagecount+1
					primers[primer]=x
				}
				if(check2){
					x = primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]
					x.phagecount=x.phagecount+1
					primers[twoBitEncoding(RevComplement(twoBitDecode(primer)))]=x
				}

			}
			//if(cluster=="A1"&&(twoBitDecode(primer)=="TGAGAGCCCCGTAGACGG")){
			//	println(primers[primer].phagecount)
			//}

		}
		//fmt.Println(phage)
	}
	primerChan<-primers
}
func makePrimers(seq string,bps int, c chan map[uint64]bool ) {
	phagprimers:=make(map[uint64]bool)
	for k:= 0; k <= len(seq) - bps; k++ {

		primer:=twoBitEncoding(seq[k:k+bps])
		phagprimers[primer]=true

	}
	c<-phagprimers
}
func DoPrimerAnalysis(from int, to int){
	f, err := os.Create(WorkingDir +"\\Data\\Unique.csv")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	for i:=from;i<=to;i++{
		//primerAnalysis(i,w)
		//primerAnalysisDouble(i,w)
		PrimerAnalysisParallel(i,w)
	}
	err = w.Flush() // Don't forget to flush!
	if err != nil {
		log.Fatal(err)
	}
}
