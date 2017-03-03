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
	sort.Ints(f)
	sort.Ints(r)

	index =0;
	//        int count =0;
	primerMatch:=make(map[PrimerMatch]float64)
	var b,frag int;
	for i:=0;i<len(f);i++{
		if(index>=len(r)){
			break
		}
		a:=f[i]
		//            System.out.fmt.Println(count);
		//            count++;
		if(index>=len(r)){
			println()
		}
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