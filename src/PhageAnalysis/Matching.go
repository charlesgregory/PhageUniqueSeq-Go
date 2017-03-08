package PhageAnalysis

import (
	"fmt"
	"math"
	"sort"
)

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
	locations := make(map[int]uint64)
	var part,part2,rprimer uint64;
	var integers,integersr []int;
	var sequence1,sequence2 string ;
	var check bool;
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
	f := make([]int,len(locations))
	index:=0
	for i,_:=range locations {
		f[index]=i
		index++
		//print(i)
		//print(" ")
	}
	//println()
	index=0
	//println()
	sort.Ints(f)
	count:=0
	//        int count =0;
	primerMatch:=make(map[PrimerMatch]float64)
	var b,a,frag int;
	for _,a=range f{
		for _,b=range f{
			frag=b-a
			if(frag<500){
				continue
			}else{
				break
			}
		}
		if(frag>2000){
			continue
		}else{
			pF := locations[a];
			pR := locations[b];
			_,check=primerTm[pF]
			_,check=primerTm[pR]
			if(math.Abs(primerTm[pF]-primerTm[pR])<5.0){
				match := PrimerMatch{pF,twoBitEncoding(RevComplement(twoBitDecode(pR)))};
				_,check:=primerMatch[match]
				if(!check){
					primerMatch[match]=float64(frag);
					count++
					//println("match")
				}else{
					delete(primerMatch,match)
					//println("delete")
				}
			}
		}
	}
	return primerMatch

}