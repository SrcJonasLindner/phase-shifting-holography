//** DM-script to test the Gun-Tilt API
// modify the <beamTiltPath>
// ReadGunTiltCoord.exe needs two command line paramters 
// 1st cmd line parameter: <integer: (10**digits of the coordinate you want to recieve)>
// 2nd cmd line parameter: <bool: 0= x-coordinate, 1= y-coordinate of the current gun tilt>
// The ReadGunTiltCoord will return the current gun tilt of the chosen coordinate as interger.
// The first digit codes the sign, the rest has to be divided by 10**digits
//

void getBeamTilt(number &x, number &y, number digits){

number factor= 10**digits;
string beamTiltPath="C:\\Users\\supervisor\\Desktop\\BeamTiltApi\\ReadGunTiltCoord.exe"
string command= beamTiltPath+ " " + factor + " 0" +"\n";

x= LaunchExternalProcess( command )

if (x > factor) {
	x=(x-factor)*-1;
}
x/=factor;
command= beamTiltPath+ " " + factor + " 1" +"\n";
y= LaunchExternalProcess( command )
if (y > factor) {
	y=(y-factor)*-1;
}
y/=factor;

}

{
 number x=0,y=0 // will hold the current beam tilt as float
 getBeamTilt(x,y,9) /
 result("tilt:" + x + " , " + y +"\n");
}