

string pathname;
string python_command="";
OpenDialog(pathname)

number fileID= OpenFileForReading(pathname)
number flag=1;

number ImageID
getpersistentnumbernote("Tilt Series:LastID", ImageID)
number PInit = 10	//initial frequency guess
TagGroup tgImg = ImageGetTagGroup(GetImageFromID(ImageID))
if(TagGroupDoesTagExist( tgImg, "PInit" )) tgImg.TagGroupGetTagAsNumber( "PInit", PInit )



string InputStringWrapper(number imageID){
//string filename= OpenDialog(pathname)
//string InputString="filename="+filename+"n"; 
string InputString="\tT="+PInit+"\n\tplot_fits=0\n\tInputID="+ imageID+ "\n"
result("\n"+inputstring)
	return InputString 
}



string SkipInput(number FileID){
	string args="", skipped;
	number skipflag=1
	string line=""
	while(skipflag==1 && line!= "\t#<InputEnd>\n"){
		skipflag= ReadFileLine(fileID,line)
		skipped+=line + "\n"
	}
	return skipped;
}


string line=""
while(flag==1){
	flag= ReadFileLine(fileID,line)
	
	if (line =="\t#<InputStart>\n"){
		result("replaced Input")
		SkipInput(fileID)
		line=InputStringWrapper(imageID);
	}
	python_command+=line;
}
closeFile(fileID)

result("\n Image ID :"+imageID)
ExecutePythonScriptString(python_command,0)
result("python done\n")