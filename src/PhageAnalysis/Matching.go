package PhageAnalysis

import (
	"fmt"
	"math"
	"sort"
	"os"
	"log"
	"bufio"
	"strconv"
)
func MatchPrimersParallel(threads int){
	f, err := os.Create(WorkingDir +"Data"+pathslash+"Matched.csv")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	phageList:=ParsePhages();
	for strain,clusters:=range phageList {

		/**
		 FOR EACH STRAIN
		 */
		fmt.Println("Starting:" + strain)
		for cluster,phages:=range clusters {

			/**
			 FOR EACH CLUSTER
			 */

			//fmt.Println(size)
			if (len(phages) > 1) {

				/**
				 GRAB PRIMERS
				 */
				fmt.Println(cluster)
				primerTm:=make(map[uint64]float64)
				matchedPrimers:=make(map[PrimerMatch][]float64)
				primers:= ReadUniquePrimers(cluster,strain)
				//var size bool=true;
				for i:=0;i<len(primers);i++{
					primerTm[primers[i]]=easytm(twoBitDecode(primers[i]))
				}
				fmt.Println(len(primers));
				fmt.Println(len(phages));
				var primChan= make(chan map[PrimerMatch]float64,len(phages))
				var jobs = make(chan MatchingWork, len(clusters))
				for w := 1; w <= threads; w++ {
					go MatchingWorker(jobs, primChan)
				}
				for _,seq:=range phages {
					jobs <- MatchingWork{seq,primers,primerTm}
				}
				close(jobs)
				count:=0;
				for range phages {
					/**
					 * FOR EACH PHAGE
					 */
					if(count==0){
						batch:=<-primChan
						for primerM,frag:=range batch{
							temp:=make([]float64,1)
							temp[0]=frag
							matchedPrimers[primerM]=temp
						}
					}else{
						batch:=<-primChan
						for primerM,frag:=range batch{
							arr,check:=matchedPrimers[primerM]
							if(check){
								temp:=make([]float64,len(arr)+1)
								for i:=0;i<len(arr);i++{
									temp[i]=arr[i]
								}
								temp[len(arr)]=frag
								matchedPrimers[primerM]=temp
							}
						}
					}
					count++
				}
				ind:=0
				for pM,arr:=range matchedPrimers{
					if(len(arr)==len(phages)){
						var sum float64 = 0.0
						var avg float64 = 0.0
						var stddev float64 = 0.0
						for _,i:=range arr{
							sum=sum+i
						}
						avg=sum/float64(len(arr))
						for _,i:=range arr{
							stddev=stddev+(math.Pow(i-avg,2))
						}
						stddev=stddev/float64(len(arr))
						stddev=math.Sqrt(stddev)
						w.WriteString(strain+",")
						w.WriteString(cluster)
						w.WriteString(",")
						w.WriteString(RevComplement(twoBitDecode(pM.F)))
						w.WriteString(",")
						w.WriteString(RevComplement(twoBitDecode(pM.R)))
						w.WriteString(",")
						w.WriteString(strconv.FormatFloat(avg, 'f',3,64))
						w.WriteString(",")
						w.WriteString(strconv.FormatFloat(stddev, 'f',3,64)+"\n")
						ind++
					}
				}
				fmt.Println("Matches Compiled");
				fmt.Println(ind);
			}
			//for primerM,arr :=range matchedPrimers{
			//	arr = matchFrags.get(m);
			//	newA = new double[arr.length];
			//	for(int i=0;i<arr.length;i++){
			//		newA[i]=arr[i];
			//	}
			//	db.insertMatchedPrimer(m.foward,m.reverse,z,x,newA);
			//	count++;
			//}
			//System.out.fmt.Println(count);
			//System.out.fmt.Println();
			//log.fmt.Println(z);
			//log.flush();
			//System.gc();
			//db.insertMatchedPrimerCommit();
		}
		//System.out.fmt.Println((System.nanoTime() - time) / Math.pow(10, 9) / 60.0);
	}
	fmt.Println("Matches Submitted");
	err = w.Flush() // Don't forget to flush!
	if err != nil {
		log.Fatal(err)
	}
}
func easytm(primer string)float64{
	a :=0;
	c :=0;
	g := 0;
	t := 0;
	var re float64;
	for i:=0;i<len(primer);i++{
		x:=primer[i]
		if(x=='A'||x=='a'){
			a++;
		}else if (x=='G'||x=='g'){
			g++;
		}else if (x=='C'||x=='c'){
			c++;
		}else if (x=='T'||x=='t'){
			t++;
		}
	}
	re=64.9 +41*(float64(g+c)-16.4)/float64(a+t+g+c)
	return re
}
func matchRoutine(seq string,primers []uint64,primerTm map[uint64]float64,
pMchan chan map[PrimerMatch]float64){
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
	pMchan<-primerMatch

}