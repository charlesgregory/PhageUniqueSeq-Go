package PhageAnalysis

import (
	"fmt"
	"time"
)
type AnalysisWork struct{
	phages map[string]string
	clusterNum uint8
	bps int
	cluster string
}
func AnalysisWorker(jobs <-chan AnalysisWork, results chan<- map[uint64]Primer) {
	for j := range jobs {
		results <- DoAnalysisWork(j.phages,j.clusterNum,j.bps,j.cluster)
	}
}
func DoAnalysisWork(phages map[string]string,
clusterNum uint8,bps int,
cluster string)map[uint64]Primer{
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
	return primers
}
func worker(id int, jobs <-chan int, results chan<- int) {
	for j := range jobs {
		fmt.Println("worker", id, "started  job", j)
		time.Sleep(time.Second)
		fmt.Println("worker", id, "finished job", j)
		results <- j * 2
	}
}
func Run() {
	jobs := make(chan int, 100)
	results := make(chan int, 100)
	for w := 1; w <= 3; w++ {
		go worker(w, jobs, results)
	}
	for j := 1; j <= 10; j++ {
		jobs <- j
	}
	close(jobs)
	for a := 1; a <= 10; a++ {
		<-results
	}
}