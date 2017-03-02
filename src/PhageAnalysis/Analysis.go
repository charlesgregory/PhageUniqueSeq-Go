package PhageAnalysis

import (
	"bufio"
	//"time"
	//"fmt"
	//"strconv"
	"math"
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
					_, check2 = primers[twoBitEncoding(revComplement(twoBitDecode(primer)))]
					if (!check&&!check2) {
						var x = make([]uint8, 1)
						x[0] = clusterNum
						primers[primer] = Primer{x, 1}
					} else {
						var x Primer;
						if (check) {
							x = primers[primer]
						} else {
							x = primers[twoBitEncoding(revComplement(twoBitDecode(primer)))]
						}
						x.phagecount=x.phagecount+1
						var found bool = false
						for i := 0; i < len(x.clusters); i++ {
							if (clustersMap[x.clusters[i]] == cluster) {
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
							primers[twoBitEncoding(revComplement(twoBitDecode(primer)))]=x
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
func MatchPrimers(){
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
			fmt.Println(cluster)
			primerTm:=make(map[uint64]float64)
			matchedPrimers:=make(map[PrimerMatch][]float64)
			primers:= ReadUniquePrimers(cluster,strain)
			//var size bool=true;
			for i:=0;i<len(primers);i++{
				primerTm[primers[i]]=easytm(twoBitDecode(primers[i]))
				//if(len(twoBitDecode(primers[i]))!=18){
				//	size=false
				//}
			}
			//fmt.Println(size)
			if (len(phages) > 1) {

				/**
				 GRAB PRIMERS
				 */

				fmt.Println(len(primers));
				fmt.Println(len(phages));
				count:=0;
				for _,seq:=range phages {
					/**
					 * FOR EACH PHAGE
					 */
					if(count==0){
						batch:=match(seq,primers,primerTm);
						for primerM,frag:=range batch{
							temp:=make([]float64,1)
							temp[0]=frag
							matchedPrimers[primerM]=temp
						}
					}else{
						batch:=match(seq,primers,primerTm);
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
			}

			fmt.Println("Matches Compiled");
			fmt.Println(len(matchedPrimers));
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
	//db.db.close();
}
func makePrimers(seq string,bps int, c chan map[uint64]bool ) {
	phagprimers:=make(map[uint64]bool)
	for k:= 0; k <= len(seq) - bps; k++ {

		primer:=twoBitEncoding(seq[k:k+bps])
		phagprimers[primer]=true

	}
	c<-phagprimers
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
func match(seq string,primers []uint64,primerTm map[uint64]float64)map[PrimerMatch]float64{
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
	forward := make(map[int]uint64)
	reverse := make(map[int]uint64)
	var part,part2,rprimer uint64;
	var integers,integersr []int;
	var sequence1,sequence2 string ;
	var check bool;
	for _,primer:=range primers{
		sequence1 = twoBitDecode(primer);
		if(len(sequence1)<10){
			fmt.Println(sequence1)
		}
		part = twoBitEncoding(sequence1[0:10]);
		integers,check= seqInd[part];
		if (check) {
			for _,num:=range integers{
				if ((len(sequence1) + num) < len(seq) &&
					seq[num:len(sequence1) + num]==sequence1) {
					forward[num]=primer;
				}
			}
		}
		rprimer = twoBitEncoding(revComplement(twoBitDecode(primer)));
		sequence2 = twoBitDecode(rprimer);
		part2 = twoBitEncoding(sequence2[0:10])
		integersr, check= seqInd[part2];
		if (check) {
			for _,num:=range integersr{
				if ((len(sequence2) + num) < len(seq) &&
					seq[num:len(sequence2) + num]==sequence2) {
					reverse[num]=rprimer;
				}
			}
		}
	}
	/**
	 * FRAGMENT FINDING
	 */
	f := make([]int,len(forward))
	r := make([]int,len(reverse))
	index:=0
	for i,_:=range forward{
		f[index]=i
		index++
	}
	index=0
	for i,_:=range reverse{
		r[index]=i
		index++
	}

	index =0;
	//        int count =0;
	primerMatch:=make(map[PrimerMatch]float64)
	var b,frag int;
	for i:=0;i<len(f);i++{
		a:=f[i]
		//            System.out.fmt.Println(count);
		//            count++;
		b=r[index];
		for(index<len(r)-1&&b<a){
			index++;
			b=r[index];
		}
		frag =b-a;
		for(frag<500&&index<len(r)-1){
			index++;
			b=r[index];
			frag = b-a;
		}
		for(frag<=2000&& index<len(r)-1){
			pF := forward[a];
			pR := reverse[b];
			if(math.Abs(primerTm[pF]-primerTm[pR])<5.0){
				match := PrimerMatch{pF,pR};
				_,check:=primerMatch[match]
				if(!check){
					primerMatch[match]=float64(frag);
				}else{
					delete(primerMatch,match)
				}
				index++;
				b=r[index];
				frag = b-a;
			}else{
				index++;
				b=r[index];
				frag = b-a;
			}
		}
		index =0;
	}
	return primerMatch

}
func DoPrimerAnalysis(){
	f, err := os.Create(workingDir +"\\Data\\Unique.csv")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	primerAnalysis(18,w)
	err = w.Flush() // Don't forget to flush!
	if err != nil {
		log.Fatal(err)
	}
}
