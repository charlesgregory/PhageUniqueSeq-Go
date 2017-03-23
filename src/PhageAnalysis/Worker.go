package PhageAnalysis

import (
	"fmt"
	"sort"
	"math"
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
type MatchingWork struct{
	seq string
	primers []uint64
	primerTm map[uint64]float64
}
func MatchingWorker(jobs <-chan MatchingWork, results chan<- map[PrimerMatch]float64) {
	for j := range jobs {
		results <- DoMatchingWork(j.seq,j.primers,j.primerTm)
	}
}
func DoMatchingWork(seq string,primers []uint64,primerTm map[uint64]float64)map[PrimerMatch]float64{
	seqInd:=make(map[uint64][]int)

	/**
	 HASH SEQ
	 */
	var rc,temp,list []int;
	var sub uint64;
	for i := 0; i <= len(seq) - 10; i++ {
		sub = twoBitEncoding(seq[i:i + 10]);
		_,check:=seqInd[sub]
		if (check) {
			rc = seqInd[sub];
			temp = make([]int,len(rc)+1);
			for j:=0;j<len(rc);j++{
				temp[j]=rc[j]
			}
			temp[len(rc)]=i;
			seqInd[sub]=temp;
		} else {
			list =make([]int,1);
			list[0]=i;
			seqInd[sub]= list;
		}
	}

	/**
	 LOCATION INDEXING
	 */
	locations := make(map[int]uint64)
	var part,part2,rprimer,pF,pR uint64;
	var integers,integersr []int;
	var sequence1,sequence2 string ;
	var check,check2 bool;
	for _,primer:=range primers{
		sequence1 = twoBitDecode(primer);
		part = twoBitEncoding(sequence1[0:10]);
		integers,check= seqInd[part];
		if (check) {
			for _,num:=range integers{
				if ((len(sequence1) + num) < len(seq) &&
					seq[num:len(sequence1) + num]==sequence1) {
					locations[num]=primer;
				}
			}
		}
		rprimer = twoBitEncoding(RevComplement(twoBitDecode(primer)));
		sequence2 = twoBitDecode(rprimer);
		part2 = twoBitEncoding(sequence2[0:10])
		integersr, check= seqInd[part2];
		if (check) {
			for _,num:=range integersr{
				if ((len(sequence2) + num) < len(seq) &&
					seq[num:len(sequence2) + num]==sequence2) {
					locations[num]=primer;
				}
			}
		}
	}
	/**
	 * FRAGMENT FINDING
	 */
	var f = make([]int,len(locations))
	index:=0
	for i,_:=range locations {
		f[index]=i
		index++
	}
	index=0
	sort.Ints(f)
	count:=0
	primerMatch:=make(map[PrimerMatch]float64)
	remove:=make(map[PrimerMatch]bool)
	var b,a,frag,i,j,llimit,ulimit int;
	var match,match2 PrimerMatch
	for i,a=range f{
		llimit=a+500
		ulimit=a+2000
		for j=i+1;j<len(f);j++{
			b=f[j]
			//println(b)
			if b>=llimit&&b<=ulimit{
				frag=b-a
				pF = locations[a];
				pR = locations[b];
				if(math.Abs(primerTm[pF]-primerTm[pR])<5.0){
					match = PrimerMatch{pF,twoBitEncoding(RevComplement(twoBitDecode(pR)))};
					_,check=primerMatch[match]
					match2 = PrimerMatch{twoBitEncoding(RevComplement(twoBitDecode(pF))),pR};
					_,check2=primerMatch[match]
					if(!check&&!check2){
						primerMatch[match]=float64(frag);
						count++
						//println("match")
					}else{
						if(check){
							remove[match]=true
						}else{
							remove[match2]=true
						}
					}
				}
			}else if b>ulimit{
				break
			}
		}
	}
	for p,_:=range remove{
		delete(primerMatch,p)
	}
	return primerMatch

}