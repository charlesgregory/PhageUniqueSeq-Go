package PhageAnalysis

import (
	"strings"
	"io/ioutil"
	"net/http"
	"fmt"
	"os"
	"io"
	"bufio"
	"strconv"
	"regexp"
)
var down=false
var multi=false
var online=false
var verbose=false
var workingDir =""

func importFile(file string)string{
	fileBytes, _:=ioutil.ReadFile(file)
	return string(fileBytes)
}
func DownloadFromUrl(url string, fileName string) {
	if(verbose){
		fmt.Println("Downloading", url, "to", fileName)
	}

	// TODO: check file existence first with io.IsExist
	output, err := os.Create(fileName)
	if err != nil {
		fmt.Println("Error while creating", fileName, "-", err)
		return
	}
	defer output.Close()

	response, err := http.Get(url)
	if err != nil {
		fmt.Println("Error while downloading", url, "-", err)
		return
	}
	defer response.Body.Close()

	n, err := io.Copy(output, response.Body)
	if err != nil {
		fmt.Println("Error while downloading", url, "-", err)
		return
	}
	if(verbose){
		fmt.Println(n, "bytes downloaded.")
	}
}
func ReadUniquePrimers(cluster string, strain string)[]uint64{
	f, _ := os.Open(workingDir +"\\Data\\Unique.csv")
	var lines []uint64
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line:=scanner.Text()
		if(strings.Contains(line,strain)&&strings.Contains(line,cluster)){
			strA:=strings.Split(line,",")
			i,_:=strconv.ParseUint(strA[2],10,64)
			lines = append(lines,i )
			//if(len(strA[2])!=18){
			//	println(strain+" "+cluster+" "+strA[2])
			//}
		}

	}
	if err := scanner.Err(); err != nil {
		fmt.Fprintln(os.Stderr, err)
	}
	return lines
}
func readUniqueClusters()map[string]map[string]bool{
	clustList:= make(map[string]map[string]bool)
	f, _ := os.Open(workingDir +"\\Data\\Unique.csv")
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line:=scanner.Text()
		strA:=strings.Split(line,",")
		_,check:=clustList[strA[0]]
		if(check){
			clustList[strA[0]][strA[1]]=true
		}else{
			clustList[strA[0]]=make(map[string]bool)
			clustList[strA[0]][strA[1]]=true
		}

	}
	if err := scanner.Err(); err != nil {
		fmt.Fprintln(os.Stderr, err)
	}
	return clustList
}

func ParseArgs(args []string)  {
	for i:=0;i<len(args);i++{
		if(args[i]=="-download") {
			down=true
		}
		if(args[i]=="-onlyMulti") {
			multi=true
		}
		if(args[i]=="-online") {
			online=true
		}
		if(args[i]=="-verbose") {
			verbose=true
		}
		if(args[i]=="-path") {
			workingDir =args[i+1]
		}
	}
	if(down){
		DownloadFromUrl("http://phagesdb.org/media/Mycobacteriophages-All.fasta", workingDir +"\\Fastas\\Mycobacteriophages-All.fasta")
		if(!multi){
			DownloadFromUrl("http://phagesdb.org/data/?set=seq&type=simple", workingDir +"\\Fastas\\PhageDBData.txt")
			DownloadFromUrl("http://phagesdb.org/data/?set=seq&type=full", workingDir +"\\Fastas\\PhageDBDataFull.txt")
		}

	}
}
func ParsePhages()map[string]map[string]map[string]string{
	phageList:=make(map[string]map[string]map[string]string)
	total:=""
	count:=1
	nextFilePat, _ := regexp.Compile("\"next\":")
	phagePat, _ := regexp.Compile("\"phage_name\":")
	clusterPat, _ := regexp.Compile("\"cluster\":")
	subclusterPatExists, _ := regexp.Compile("\"psubcluster\":")
	subclusterPat, _ := regexp.Compile("\"subcluster\":")
	strainPat, _ := regexp.Compile("\"genus\":")
	fastaPat, _ := regexp.Compile("\"fasta_file\":")
	var nextFile="http://phagesdb.org/api/sequenced_phages/"
	var firstLine string
	var scanner *bufio.Reader
	var file *os.File
	var line []byte
	var name,strain,cluster,subclusterExists, subcluster,fasta string
	var check bool
	for nextFile!="ul"{
		if(online) {
			DownloadFromUrl(nextFile,
				workingDir +"\\Fastas\\Phagelist" + strconv.Itoa(count) + ".txt")
		}
		file,_=os.Open(workingDir +"\\Fastas\\Phagelist"+strconv.Itoa(count)+".txt")
		scanner = bufio.NewReader(file)
		line,_,_=scanner.ReadLine()
		firstLine=string(line)
		for _,x:=range nextFilePat.FindAllStringIndex(firstLine,-1){
			sub:=firstLine[x[0]+(x[1]-x[0])+1:]
			nextFile=sub[:strings.Index(sub,",")-1]
		}
		var str string=firstLine
		for len(line)>0{
			line,_,_=scanner.ReadLine()
			str=str+string(line)
		}
		total=total+str
		count++
	}
	for _,x:=range phagePat.FindAllStringIndex(total,-1){
		sub:=total[x[0]+(x[1]-x[0])+1:]
		name=sub[:strings.Index(sub,",")-1]
		x:= clusterPat.FindStringIndex(sub)
		sub=sub[x[0]+(x[1]-x[0])+1:]
		cluster =sub[:strings.Index(sub,",")-2]
		x= subclusterPatExists.FindStringIndex(sub)
		sub=sub[x[0]+(x[1]-x[0]):]
		subclusterExists =sub[:strings.Index(sub,",")]
		if(subclusterExists=="null"){
			subcluster=cluster
		}else{
			x= subclusterPat.FindStringIndex(sub)
			sub=sub[x[0]+(x[1]-x[0])+1:]
			subcluster =sub[:strings.Index(sub,",")-2]
		}
		x=strainPat.FindStringIndex(sub)
		sub=sub[x[0]+(x[1]-x[0])+1:]
		strain=sub[:strings.Index(sub,",")-1]
		x=fastaPat.FindStringIndex(sub)
		sub=sub[x[0]+(x[1]-x[0])+1:]
		fasta=sub[:strings.Index(sub,",")-1]
		_,check=phageList[strain]
		if(!check){
			phageList[strain]=make(map[string]map[string]string)
		}
		_,check=phageList[strain][subcluster]
		if(!check){
			phageList[strain][subcluster]=make(map[string]string)
		}
		if(subcluster=="Singleton"){
			phageList[strain][name]=make(map[string]string)
			phageList[strain][name][name]=fasta
		}else{
			phageList[strain][subcluster][name]=fasta
		}

	}
	var seq string
	for strain,clusters:=range phageList{
		for cluster,phages:=range clusters {
			for phage,fasta:=range phages {
				//println(strain+" "+cluster+" "+phage+" "+fasta)
				if(online) {
					DownloadFromUrl(fasta,
						workingDir +"\\Fastas\\" + phage + ".fasta")

				}
				seq=importFile(workingDir +"\\Fastas\\" + phage + ".fasta")
				lines:=strings.Split(seq,"\n")
				seq:=""
				for j:=1;j<len(lines) ; j++ {
					seq=seq+lines[j]
				}
				phageList[strain][cluster][phage]=seq

			}
		}
	}
	return phageList
}
