package PhageAnalysis

type Primer struct{
	clusters map[uint8]uint8
}
func (p *Primer)addCluster(c uint8){
	var temp map[uint8]uint8
	if(len(p.clusters)>0){
		p.clusters[c]=0
	}else{
		temp=make(map[uint8]uint8)
		temp[c]=0
		p.clusters=temp
	}
}
func (p *Primer)addPhage(c uint8){
	p.clusters[c]=p.clusters[c]+1
}
func (p *Primer)combine(prim Primer){
	var key uint8
	var value uint8
	var check bool
	for key,value=range prim.clusters{
		_,check = p.clusters[key]
		if(check){
			p.clusters[key]=p.clusters[key]+value
		}else{
			p.clusters[key]=value
		}

	}
}
type PrimerMatch struct {
	F uint64
	R uint64
}