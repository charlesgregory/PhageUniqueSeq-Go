package PhageAnalysis
var ENCODING_MASK uint64 = 3
var LEN_MASK uint64=31
func twoBitEncoding(seq string )uint64{
	var encoded uint64= 0
	var curShift uint64=5
	encoded =encoded ^ uint64(len(seq))
	for i:=0;i<len(seq);i++{
		var baseEncoded uint64
		if(seq[i]=='A'){
			baseEncoded=0
		}else if (seq[i]=='G'){
			baseEncoded=1
		}else if (seq[i]=='C'){
			baseEncoded=2
		}else {
			baseEncoded=3
		}
		encoded=baseEncoded << curShift ^ encoded
		curShift+=2
	}
	return encoded
}
func twoBitDecode(encoded uint64)string{
	var length uint64= encoded & LEN_MASK
	builder :=""
	var i uint64=5;
	for ;i<=62;i=i+2{
		if((i-5)/2==length){
			break;
		}
		var maskedLetter uint64= encoded & (ENCODING_MASK << i);
		var decodedChar uint64 = maskedLetter >> i;
		var newChar string;
		if(decodedChar==0){
			newChar="A";
		}else if(decodedChar==1){
			newChar="G";
		}else if(decodedChar==2){
			newChar="C";
		}else{
			newChar="T";
		}
		builder=builder+newChar;

	}
	return builder;
}
func revComplement(seq string)string{
	newString:=""
	for i:=0;i<len(seq);i++{
		if(seq[len(seq)-1-i]=='G'||seq[len(seq)-1-i]=='g'){
			newString=newString+"C"
		}else if(seq[len(seq)-1-i]=='C'||seq[len(seq)-1-i]=='c'){
			newString=newString+"G"
		}else if(seq[len(seq)-1-i]=='A'||seq[len(seq)-1-i]=='a'){
			newString=newString+"T"
		}else if(seq[len(seq)-1-i]=='T'||seq[len(seq)-1-i]=='t'){
			newString=newString+"A"
		}
	}
	return newString
}
func reEncodeReverseComplementTwoBit(encoding uint64 )uint64 {
	var length uint64= encoding & LEN_MASK;
	var newEncoding uint64= 0;
	newEncoding = newEncoding | length;
	var i uint64=5;
	for ;i<=62;i=i+2 {
		if((i-5)/2==length){
			break;
		}
		//Only have the bits of the letter
		var maskedLetter uint64= encoding & (ENCODING_MASK << i);
		var rcEncoded uint64= maskedLetter >> i;
		if(rcEncoded==0) {
			rcEncoded = 3;
		}else if (rcEncoded==1) {
			rcEncoded = 2;
		}else if (rcEncoded==2) {
			rcEncoded = 1;
		}else{
			rcEncoded=0;
		}
		newEncoding = (rcEncoded << ((length)+4)*2-i | newEncoding);
	}

	return newEncoding;
}
