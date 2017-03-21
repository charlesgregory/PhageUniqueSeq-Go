package PhageAnalysis

type Primer struct{
	clusters map[uint8]bool
	phagecount uint16
}
func (p *Primer)addCluster(c uint8){
	var temp map[uint8]bool
	if(len(p.clusters)>0){
		p.clusters[c]=true
	}else{
		temp=make(map[uint8]bool)
		temp[c]=true
		p.clusters=temp
	}
}
func (p *Primer)addPhage(){
	p.phagecount++
}
func (p *Primer)combine(prim Primer){
	var key uint8
	for key,_=range prim.clusters{
		p.clusters[key]=true
	}
	p.phagecount=p.phagecount+prim.phagecount
}
type PrimerMatch struct {
	F uint64
	R uint64
}