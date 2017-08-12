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
	//"time"
	"math"
	"sort"
	"regexp"
	"sync"
	"log"
	"github.com/texttheater/golang-levenshtein/levenshtein"
	"time"
)

func readTest()[]uint64{
	f, _ := os.Open(WorkingDir +""+pathslash+"Data"+pathslash+"Test.csv")
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
	}
	count:=0
	for _,seq:=range phages{
		//for _,p:=range primers{
		//	ind:=strings.Index(seq,twoBitDecode(p))
		//	if(ind==-1){
		//		ind=strings.Index(RevComplement(seq),RevComplement(twoBitDecode(p)))
		//		if(ind==-1){
		//			count++
		//		}
		//	}
		//}
		//println(count)
		if(count==0){
			batch:=matchTest(seq,primers,primerTm)
			for primerM,frag:=range batch{
				fmt.Print(twoBitDecode(primerM.F)+" "+twoBitDecode(primerM.R))
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
		primers:=ReadUniquePrimers(cluster,strain)
		//println(len(primers))
		//println(len(phages))
		var re=make(chan string, len(phages))
		for phage,seq:=range phages{
			wg.Add(1)
			go testUniquePresent(primers,seq,re,phage)
			//for _,p:=range primers{
			//	if(!strings.Contains(seq,twoBitDecode(p))&&!
			//		strings.Contains(seq, RevComplement(twoBitDecode(p)))){
			//		println(phage+" "+twoBitDecode(p)+" "+RevComplement(twoBitDecode(p)))
			//	}
			//}
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
				println("\t"+otherc)
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
func UniqueConfirmCluster(strain string,cluster string){
	phagList:=ParsePhages()
	phages:=phagList[strain][cluster]
	println("Checking:"+cluster)
	primers:=ReadUniquePrimers(cluster,strain)
	println(len(primers))
	println(len(phages))
	var re=make(chan string, len(phages))
	for phage,seq:=range phages{
		wg.Add(1)
		go testUniquePresent(primers,seq,re,phage)
		//for _,p:=range primers{
		//	if(!strings.Contains(seq,twoBitDecode(p))&&!
		//		strings.Contains(seq, RevComplement(twoBitDecode(p)))){
		//		println(phage+" "+twoBitDecode(p)+" "+RevComplement(twoBitDecode(p)))
		//	}
		//}
	}
	wg.Wait()
	for range phages{
		x:=<-re
		print(x)
	}
	//println()
	//println("checking others")
	//for otherc,phagesc:=range clusters{
	//	if(otherc!=cluster){
	//		println("\t"+otherc)
	//		var re=make(chan string, len(phagesc))
	//		for _,seq:=range phagesc{
	//			wg.Add(1)
	//			go testUniqueNotPresent(primers,seq,re)
	//		}
	//		wg.Wait()
	//		for range phagesc{
	//			x:=<-re
	//			print(x)
	//		}
	//	}
	//
	//}

}
func FindPrimer(strain string,p string){
	phagList:=ParsePhages()
	clusters:=phagList[strain]
	phageCount:=0
	clusterMap:=make(map[string]bool)
	for cluster,phages:=range clusters{
		temp:=phageCount
		for _,seq:=range phages{
			if(strings.Contains(seq,p)||
				strings.Contains(seq, RevComplement(p))){
				phageCount++
			}
		}
		if phageCount!=temp{
			clusterMap[cluster]=true
		}
	}
	var clust string
	for k,_:=range clusterMap{
		println(k)
		clust=k
	}
	println(len(clusters[clust]))
	println(phageCount)

}
func testUniquePresent(primers []uint64,seq string,re chan string,phage string){
	defer wg.Done()
	var found bool=true
	var str string=""
	for _,p:=range primers{
		if(!strings.Contains(seq,twoBitDecode(p))&&!
			strings.Contains(seq, RevComplement(twoBitDecode(p)))){
			str=str+twoBitDecode(p)+" "+RevComplement(twoBitDecode(p))+" "+phage+"\n"
			found=false
		}
	}
	if(found){
		re<-""
	}else{
		re<-str
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
		re<-"found"
	}else{
		re<-""
	}
}
func MatchingConfirm(strain string){
	phagList:=ParsePhages()
	clusters:=phagList[strain]
	var prim map[PrimerMatch]Stat
	for cluster,phages:=range clusters{
		prim=ReadMatchedPrimers(cluster,strain)
		var arr =make([]float64,len(phages))
		println(cluster)
		for p,s:=range prim{
			var sum float64 = 0.0
			var avg float64 = 0.0
			var stddev float64 = 0.0
			count:=0
			for _,seq:=range phages{
				locf:=strings.Index(seq,twoBitDecode(p.F))
				locr:=strings.Index(seq,RevComplement(twoBitDecode(p.R)))
				if(locf==-1||locr==-1){
					locf=strings.Index(seq,RevComplement(twoBitDecode(p.F)))
					locr=strings.Index(seq,twoBitDecode(p.R))
				}
				arr[count]=math.Abs(float64(locf)-float64(locr))
				count++
			}
			for _,i:=range arr{
				sum=sum+i
			}
			avg=sum/float64(len(arr))
			for _,i:=range arr{
				stddev=stddev+(math.Pow(i-avg,2))
			}
			stddev=stddev/float64(len(arr))
			stddev=math.Sqrt(stddev)
			if(math.Abs(s.Mean-avg)>0.01||math.Abs(s.Stddev-stddev)>0.01){
				fmt.Print(cluster)
				fmt.Print(" ")
				fmt.Print(twoBitDecode(p.F))
				fmt.Print(" ")
				fmt.Print(twoBitDecode(p.R))
				fmt.Print(" ")
				fmt.Print(s.Mean)
				fmt.Print(" ")
				fmt.Print(avg)
				fmt.Print(" ")
				fmt.Print(s.Stddev)
				fmt.Print(" ")
				fmt.Print(stddev)
				fmt.Println()
			}
		}
	}
}
func testMatchPresent(){}
func TestRegex(){
	seq:=ReadFile(WorkingDir+"Fastas"+pathslash+"20ES.fasta")
	x,_:=regexp.Compile("GATCGTC")
	for _,y:=range x.FindAllStringIndex(seq,-1){
		print(y[0])
		print(" ")
		println(y[1])
	}

}
func TestPrimer(){
	var x Primer
	x.addCluster(1)
	println(len(x.clusters))

}
func ExportClusterSummary(){
	f, err := os.Create(WorkingDir +"Data"+pathslash+"ClusterSummary.csv")
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	phageList:=ParsePhages()
	for strain,clusters:=range phageList{
		for cluster,phages:=range clusters{
			primers:=ReadUniquePrimers(cluster,strain)
			w.WriteString(strain+",")
			w.WriteString(cluster)
			w.WriteString(",")
			w.WriteString(strconv.Itoa(len(phages)))
			w.WriteString(",")
			w.WriteString(strconv.Itoa(len(primers))+"\n")
		}
	}
	err = w.Flush() // Don't forget to flush!
	if err != nil {
		log.Fatal(err)
	}
}
func Test(){
	t1:=time.Now()
	for i:=0;i<10000000;i++{
		levenshtein.DistanceForStrings(
			[]rune("GATGATGCTAGCGATCGA"),[]rune("GATTGATCGGACTAGGCT"),
			levenshtein.DefaultOptions)
	}
	fmt.Println(time.Since(t1).Seconds())
	//for z,_:=range y{
	//	println(twoBitDecode(z))
	//}
}