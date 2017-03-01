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
)
var down=false
var multi=false
var online=false

type Fasta struct{
	header string
	seq string
}
func getPhageOnlineinfo(header []string)[2]string{
	name:=header[2]
	var phage [2]string
	url:="http://phagesdb.org/phages/"+name+"/"
	response, err := http.Get(url)
	if err != nil {
		fmt.Println("Error while downloading", url, "-", err)
	}
	defer response.Body.Close()
	bodyBytes, _ := ioutil.ReadAll(response.Body)
	res:=string(bodyBytes)
	strainStr:="<tr><td class=\"detailLabel\">Isolation Host</td><td class=\"detailValue\">"
	clustStr:="<tr><td class=\"detailLabel\">Subcluster</td><td class=\"detailValue\">"
	strain:=strings.Index(res,strainStr)+len(strainStr)
	clust:=strings.Index(res,clustStr)+len(clustStr)
	if(strain!=-1&&clust!=-1){
		finalStrain:=strings.Split(res[strain:strain+100],">")
		if(len(finalStrain)>2){
			finalStrain2:=strings.Split(finalStrain[2],"<")[0]
			finalClust:=strings.Split(strings.Split(res[clust:clust+100],">")[1],"<")[0]
			phage[0]=finalStrain2
			if(finalClust!=""){
				phage[1]=finalClust
			}else {
				phage[1]=name
			}
		}
	}
	return phage
}
func getPhageInfo(ffile string)[][3]string{
	var phages [][3]string
	var clustStrain [3]string
	lines:=strings.Split(ffile,"\n")
	for i := 0; i < len(lines); i++ {
		col:=strings.Split(lines[i],"\t")
		if(len(col)>2){
			clustStrain[0]=col[1]
			clustStrain[1]=col[2]
			clustStrain[2]=col[0]
			phages=append(phages,clustStrain)
		}

	}
	return phages
}
func importMulti() []Fasta  {
	file := importFile("C:\\Users\\musta_000\\IdeaProjects\\GoProjects\\Fastas\\Mycobacteriophages-All.fasta")
	seqs:=strings.Split(file,">")
	fastas :=make([]Fasta,len(seqs)-1,len(seqs)-1)
	for i:=1;i<len(seqs);i++  {
		lines:=strings.Split(seqs[i],"\n")
		seq:=""
		for j:=1;j<len(lines) ; j++ {
			seq=seq+lines[j]
		}
		fastas[i-1]=Fasta{lines[0],seq}
	}
	return fastas
}
func importFile(file string)string{
	fileBytes, _:=ioutil.ReadFile(file)
	return string(fileBytes)
}
func DownloadFromUrl(url string, fileName string) {
	fmt.Println("Downloading", url, "to", fileName)

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

	fmt.Println(n, "bytes downloaded.")
}
func GetFastas()map[string]map[string]map[string]string{
	phaglist:= make(map[string]map[string]map[string]string)
	ffile:=importFile("C:\\Users\\musta_000\\IdeaProjects\\GoProjects\\Fastas\\PhageDBDataFull.txt")
	rawlist:=getPhageInfo(ffile)
	newFastaList:=importMulti()
	var fasta Fasta
	var x, name,strain string
	var words[]string
	var z [3]string
	var check bool
	for _,fasta=range newFastaList{
		x=fasta.header
		if(len(words)>=3) {
			words = strings.Split(x, " ");
			//            System.out.println(x);
			strain = strings.Replace(strings.Replace(words[0], ">", "", -1), "Mycobacteriophage", "Mycobacterium", -1)
			_, check = phaglist[strain]
			if (!check) {
				phaglist[strain] = make(map[string]map[string]string)
			}
			name = strings.Replace(strings.Replace(words[2], ",", "", -1), "_complete", "", -1);
			if (name == "sequence") {
				name = strings.Replace(words[0], ",", "", -1);
			} else if (name == "Revised") {
				name = strings.Replace(words[3], ",", "", -1);
			}
			for _, z = range rawlist {
				if (strings.ToLower(z[2]) == strings.ToLower(name)) {
					_, check = phaglist[strain][z[1]]
					if (!check) {
						phaglist[strain][z[1]] = make(map[string]string)
					}
					phaglist[strain][z[1]][name] = fasta.seq
				}
			}
		}

	}
	return phaglist
}
func readUniquePrimers(cluster string, strain string)[]uint64{
	f, _ := os.Open("C:\\Users\\musta_000\\IdeaProjects\\GoProjects\\Data\\Unique.csv")
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
	f, _ := os.Open("C:\\Users\\musta_000\\IdeaProjects\\GoProjects\\Data\\Unique.csv")
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
//func createClusterPrimerMap()map[string]map[string][]uint64{
//	phaglist:= make(map[string]map[string][]uint64)
//	fastas:=importMulti()
//	for i:=0;i<len(fastas) ;i++  {
//		if(fastas[i].strain!=""&&fastas[i].cluster!=""){
//			phaglist[fastas[i].strain]=make(map[string][]uint64)
//		}
//	}
//	for i:=0;i<len(fastas) ;i++  {
//		if(fastas[i].strain!=""&&fastas[i].cluster!=""){
//			phaglist[fastas[i].strain][fastas[i].cluster]=make([]uint64,1)
//		}
//	}
//	for strain,clusters:=range phaglist{
//		for cluster,_:=range clusters{
//			phaglist[strain][cluster]= readUniquePrimers(cluster,strain)
//		}
//	}
//	return phaglist
//}
//func createPhageMap()map[string]map[string]map[string]string{
//	phaglist:= make(map[string]map[string]map[string]string)
//	fastas:=importMulti()
//	for i:=0;i<len(fastas) ;i++  {
//		if(fastas[i].strain!=""&&fastas[i].cluster!=""){
//			phaglist[fastas[i].strain]=make(map[string]map[string]string)
//		}
//	}
//	for i:=0;i<len(fastas) ;i++  {
//		if(fastas[i].strain!=""&&fastas[i].cluster!=""){
//			phaglist[fastas[i].strain][fastas[i].cluster]=make(map[string]string)
//		}
//	}
//	for i:=0;i<len(fastas) ;i++  {
//		if(fastas[i].strain!=""&&fastas[i].cluster!=""){
//			phaglist[fastas[i].strain][fastas[i].cluster]=make(map[string]string)
//		}
//	}
//	count:=0
//	for i:=0;i<len(fastas) ;i++  {
//		if(fastas[i].strain!=""&&fastas[i].cluster!=""){
//			phaglist[fastas[i].strain][fastas[i].cluster][fastas[i].name]=fastas[i].seq
//		}else{
//			count++
//		}
//	}
//	print(count)
//	println(" phages not used")
//	return phaglist
//}


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
	}
	if(down){
		DownloadFromUrl("http://phagesdb.org/media/Mycobacteriophages-All.fasta","C:\\Users\\musta_000\\IdeaProjects\\GoProjects\\Fastas\\Mycobacteriophages-All.fasta")
		if(!multi){
			DownloadFromUrl("http://phagesdb.org/data/?set=seq&type=simple","C:\\Users\\musta_000\\IdeaProjects\\GoProjects\\Fastas\\PhageDBData.txt")
			DownloadFromUrl("http://phagesdb.org/data/?set=seq&type=full","C:\\Users\\musta_000\\IdeaProjects\\GoProjects\\Fastas\\PhageDBDataFull.txt")
		}

	}
}
