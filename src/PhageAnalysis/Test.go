package PhageAnalysis

import (
	"fmt"
	"os"
	//"log"
	"bufio"
	//"time"
	"strconv"
	"strings"
	//"time"
	//"log"
	"time"
	"math"
	"sort"
	"regexp"
	"sync"
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

func PrimerAnalysisTest(bps int,strain string){
	phageList:=ParsePhages()
	//strain:="Mycobacterium"
	//clusters:=phageList[strain]
	var cluster,seq string
	var clusters map[string]map[string]string
	var phages map[string]string
	var primers map[uint64]Primer
	var phagprimers map[uint64]bool
	var clustersMap map[uint8]string
	//var primChan chan map[uint64]bool
	var primer uint64
	var check,check2 bool
	t1:=time.Now()
	primers= make(map[uint64]Primer)
	clustersMap=make(map[uint8]string)
	fmt.Println(strain)
	var clusterNum uint8=0
	clusters=phageList[strain]
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
	for _,v:=range primers{

		if(len(v.clusters)==1){
			//count++
			primerClust:=v.clusters[0]
			if(len(clusters[clustersMap[primerClust]])==v.phagecount){
				if(clustersMap[primerClust]=="A1"){
					count++
				}
			}
		}
	}
	fmt.Println(time.Since(t2).Minutes())
	fmt.Println(count)

}
func testPrimerMatchTest(){
	m:=make(map[PrimerMatch]bool)
	x:=PrimerMatch{twoBitEncoding("GTACGA"),twoBitEncoding("GATCAG")}
	y:=PrimerMatch{twoBitEncoding("GTACGA"),twoBitEncoding("GATCAG")}
	println(x==y)
	m[x]=true
	println(m[y]==true)

}
func readTest()[]uint64{
	f, _ := os.Open(WorkingDir +"\\Data\\Test.csv")
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
func newMatch(seq string,primers []uint64,primerTm map[uint64]float64)map[PrimerMatch]float64{
	primerM:=make(map[PrimerMatch]float64)
	forward:=make(map[int]uint64)
	reverse:=make(map[int]uint64)
	forwardFind,_:=regexp.Compile("")
	reverseFind,_:=regexp.Compile("")
	for _,p:=range primers{
		forwardFind,_=regexp.Compile(twoBitDecode(p))
		reverseFind,_=regexp.Compile(RevComplement(twoBitDecode(p)))
		fF:=forwardFind.FindAllStringIndex(seq,-1)
		rF:=reverseFind.FindAllStringIndex(seq,-1)
		for _,loc:=range fF{
			forward[loc[0]]=p
			//println(loc[0])
		}
		for _,loc:=range rF{
			reverse[loc[0]]=twoBitEncoding(RevComplement(twoBitDecode(p)))
			//println(loc[0])
		}

	}
	ind:=0
	var f=make([]int,len(forward))
	var r=make([]int,len(reverse))
	for i,_:=range forward{
		f[ind]=i
		ind++
	}
	ind=0
	for i,_:=range reverse{
		r[ind]=i
		ind++
	}
	sort.Ints(f)
	sort.Ints(r)
	println(len(f))
	println(len(r))
	var frag,b,a,i,j int
	for i=0;i<len(f);i++{
		a=f[i]
		for j=0;j<len(r);j++{
			b=r[j]
			frag=b-a
			if(!(frag<500)){
				break
			}
		}
		if(frag>2000){
			continue
		}else{
			if (math.Abs(primerTm[forward[a]]-primerTm[reverse[b]])<5.0) {
				match:=PrimerMatch{forward[a],reverse[b]}
				primerM[match]=float64(frag)
			}
		}
	}
	return primerM
}
func matchTest(seq string,primers []uint64,primerTm map[uint64]float64)map[PrimerMatch]float64{
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
		rprimer = twoBitEncoding(RevComplement(twoBitDecode(primer)));
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
		//print(i)
		//print(" ")
	}
	//println()
	index=0
	for i,_:=range reverse{
		r[index]=i
		index++
		//print(i)
		//print(" ")
	}
	//println()
	sort.Ints(f)
	sort.Ints(r)
	count:=0
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
			_,check=primerTm[pF]
			_,check=primerTm[pR]
			if(math.Abs(primerTm[pF]-primerTm[pR])<5.0){
				match := PrimerMatch{pF,pR};
				_,check:=primerMatch[match]
				if(!check){
					primerMatch[match]=float64(frag);
					count++
					//println("match")
				}else{
					delete(primerMatch,match)
					//println("delete")
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
	println(count)
	return primerMatch

}
func MatchingTest(strain string, cluster string){
	/**
			 FOR EACH CLUSTER
			 */
	phageList:=ParsePhages()
	primerTm:=make(map[uint64]float64)
	//matchedPrimers:=make(map[PrimerMatch][]float64)
	primers:= ReadUniquePrimers(cluster,strain)
	//var size bool=true;
	phages:=phageList[strain][cluster]
	for _,primer:=range primers{
		primerTm[primer]=easytm(twoBitDecode(primer))
		primerTm[twoBitEncoding(RevComplement(twoBitDecode(primer)))]=
			easytm(RevComplement(twoBitDecode(primer)))
	}
	count:=0
	for _,seq:=range phages{
		if(count==0){
			batch:=newMatch(seq,primers,primerTm)
			for primerM,frag:=range batch{
				fmt.Print(twoBitDecode(primerM.f)+" "+twoBitDecode(primerM.r))
				fmt.Print(" ")
				fmt.Println(frag)
				//temp:=make([]float64,1)
				//temp[0]=frag

			}
		}
	}
}
var wg sync.WaitGroup
func UniqueConfirm(strain string){
	phagList:=ParsePhages()
	clusters:=phagList[strain]
	for cluster,phages:=range clusters{
		println("Checking:"+cluster)
		primers:=ReadUniquePrimers(strain,cluster)
		var re=make(chan string, len(phages))
		for _,seq:=range phages{
			wg.Add(1)
			go testUniquePresent(primers,seq,re)
		}
		wg.Wait()
		for range phages{
			x:=<-re
			print(x)
		}
		println()
		println("checking others")
		for otherc,phagesc:=range clusters{
			if(otherc!=cluster){
				var re=make(chan string, len(phagesc))
				for _,seq:=range phagesc{
					wg.Add(1)
					go testUniqueNotPresent(primers,seq,re)
				}
				wg.Wait()
				for range phagesc{
					x:=<-re
					print(x)
				}
			}

		}
	}

}
func testUniquePresent(primers []uint64,seq string,re chan string){
	defer wg.Done()
	var found bool=true
	for _,p:=range primers{
		if(!strings.Contains(seq,twoBitDecode(p))&&!
			strings.Contains(seq, RevComplement(twoBitDecode(p)))){
			found=false
		}
	}
	if(found){
		re<-""
	}else{
		re<-"not found"
	}
}
func testUniqueNotPresent(primers []uint64,seq string,re chan string){
	defer wg.Done()
	var found bool=false
	for _,p:=range primers{
		if(strings.Contains(seq,twoBitDecode(p))||
			strings.Contains(seq, RevComplement(twoBitDecode(p)))){
			found=true
		}
	}
	if(found){
		re<-"not found"
	}else{
		re<-""
	}
}
func TestRegex(){
	seq:=ReadFile(WorkingDir+"\\Fastas\\20ES.fasta")
	x,_:=regexp.Compile("GATCGTC")
	for _,y:=range x.FindAllStringIndex(seq,-1){
		print(y[0])
		print(" ")
		println(y[1])
	}

}