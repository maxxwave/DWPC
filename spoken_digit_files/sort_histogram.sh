#!/bin/bash

cat $1 | awk 'BEGIN{c0=0; c1=0;c2=0;c3=0; c4=0; c5=0; c6=0; c7=0; c8=0; c9=0}
		{ if($2==0){    if($1==0){ c0 +=1 }
				if($1==1){ c1 +=1 }
				if($1==2){ c2 +=1 }
				if($1==3){ c3 +=1 }
				if($1==4){ c4 +=1 }
				if($1==5){ c5 +=1 }
				if($1==6){ c6 +=1 }
				if($1==7){ c7 +=1 }
				if($1==8){ c8 +=1 }
				if($1==9){ c9 +=1 }
				       	print( c0 "\n" c1 "\n" c2 "\n" c3 "\n" c4 "\n" c5 "\n" c6 "\n" c7 "\n" c8 "\n" c9 )}  }' >digit0


cat $1 | awk 'BEGIN{c0=0; c1=0;c2=0;c3=0; c4=0; c5=0; c6=0; c7=0; c8=0; c9=0}
		{ if($2==1){    if($1==0){ c0 +=1 }
				if($1==1){ c1 +=1 }
				if($1==2){ c2 +=1 }
				if($1==3){ c3 +=1 }
				if($1==4){ c4 +=1 }
				if($1==5){ c5 +=1 }
				if($1==6){ c6 +=1 }
				if($1==7){ c7 +=1 }
				if($1==8){ c8 +=1 }
				if($1==9){ c9 +=1 }
				       	print( c0 "\n" c1 "\n" c2 "\n" c3 "\n" c4 "\n" c5 "\n" c6 "\n" c7 "\n" c8 "\n" c9 )}  }' >digit1


cat $1 | awk 'BEGIN{c0=0; c1=0;c2=0;c3=0; c4=0; c5=0; c6=0; c7=0; c8=0; c9=0}
		{ if($2==2){    if($1==0){ c0 +=1 }
				if($1==1){ c1 +=1 }
				if($1==2){ c2 +=1 }
				if($1==3){ c3 +=1 }
				if($1==4){ c4 +=1 }
				if($1==5){ c5 +=1 }
				if($1==6){ c6 +=1 }
				if($1==7){ c7 +=1 }
				if($1==8){ c8 +=1 }
				if($1==9){ c9 +=1 }
				       print( c0 "\n" c1 "\n" c2 "\n" c3 "\n" c4 "\n" c5 "\n" c6 "\n" c7 "\n" c8 "\n" c9 )}  }' >digit2



cat $1 | awk 'BEGIN{c0=0; c1=0;c2=0;c3=0; c4=0; c5=0; c6=0; c7=0; c8=0; c9=0}
		{ if($2==3){    if($1==0){ c0 +=1 }
				if($1==1){ c1 +=1 }
				if($1==2){ c2 +=1 }
				if($1==3){ c3 +=1 }
				if($1==4){ c4 +=1 }
				if($1==5){ c5 +=1 }
				if($1==6){ c6 +=1 }
				if($1==7){ c7 +=1 }
				if($1==8){ c8 +=1 }
				if($1==9){ c9 +=1 }
				       print( c0 "\n" c1 "\n" c2 "\n" c3 "\n" c4 "\n" c5 "\n" c6 "\n" c7 "\n" c8 "\n" c9 )}  }' >digit3



cat $1 | awk 'BEGIN{c0=0; c1=0;c2=0;c3=0; c4=0; c5=0; c6=0; c7=0; c8=0; c9=0}
		{ if($2==4){    if($1==0){ c0 +=1 }
				if($1==1){ c1 +=1 }
				if($1==2){ c2 +=1 }
				if($1==3){ c3 +=1 }
				if($1==4){ c4 +=1 }
				if($1==5){ c5 +=1 }
				if($1==6){ c6 +=1 }
				if($1==7){ c7 +=1 }
				if($1==8){ c8 +=1 }
				if($1==9){ c9 +=1 }
				       	print( c0 "\n" c1 "\n" c2 "\n" c3 "\n" c4 "\n" c5 "\n" c6 "\n" c7 "\n" c8 "\n" c9 )}  }' >digit4


cat $1 | awk 'BEGIN{c0=0; c1=0;c2=0;c3=0; c4=0; c5=0; c6=0; c7=0; c8=0; c9=0}
		{ if($2==5){    if($1==0){ c0 +=1 }
				if($1==1){ c1 +=1 }
				if($1==2){ c2 +=1 }
				if($1==3){ c3 +=1 }
				if($1==4){ c4 +=1 }
				if($1==5){ c5 +=1 }
				if($1==6){ c6 +=1 }
				if($1==7){ c7 +=1 }
				if($1==8){ c8 +=1 }
				if($1==9){ c9 +=1 }
				       	print( c0 "\n" c1 "\n" c2 "\n" c3 "\n" c4 "\n" c5 "\n" c6 "\n" c7 "\n" c8 "\n" c9 )}  }' >digit5


cat $1 | awk 'BEGIN{c0=0; c1=0;c2=0;c3=0; c4=0; c5=0; c6=0; c7=0; c8=0; c9=0}
		{ if($2==6){    if($1==0){ c0 +=1 }
				if($1==1){ c1 +=1 }
				if($1==2){ c2 +=1 }
				if($1==3){ c3 +=1 }
				if($1==4){ c4 +=1 }
				if($1==5){ c5 +=1 }
				if($1==6){ c6 +=1 }
				if($1==7){ c7 +=1 }
				if($1==8){ c8 +=1 }
				if($1==9){ c9 +=1 }
				       	print( c0 "\n" c1 "\n" c2 "\n" c3 "\n" c4 "\n" c5 "\n" c6 "\n" c7 "\n" c8 "\n" c9 )}  }' >digit6



cat $1 | awk 'BEGIN{c0=0; c1=0;c2=0;c3=0; c4=0; c5=0; c6=0; c7=0; c8=0; c9=0}
		{ if($2==7){    if($1==0){ c0 +=1 }
				if($1==1){ c1 +=1 }
				if($1==2){ c2 +=1 }
				if($1==3){ c3 +=1 }
				if($1==4){ c4 +=1 }
				if($1==5){ c5 +=1 }
				if($1==6){ c6 +=1 }
				if($1==7){ c7 +=1 }
				if($1==8){ c8 +=1 }
				if($1==9){ c9 +=1 }
				       	print( c0 "\n" c1 "\n" c2 "\n" c3 "\n" c4 "\n" c5 "\n" c6 "\n" c7 "\n" c8 "\n" c9 )}  }' >digit7


cat $1 | awk 'BEGIN{c0=0; c1=0;c2=0;c3=0; c4=0; c5=0; c6=0; c7=0; c8=0; c9=0}
		{ if($2==8){    if($1==0){ c0 +=1 }
				if($1==1){ c1 +=1 }
				if($1==2){ c2 +=1 }
				if($1==3){ c3 +=1 }
				if($1==4){ c4 +=1 }
				if($1==5){ c5 +=1 }
				if($1==6){ c6 +=1 }
				if($1==7){ c7 +=1 }
				if($1==8){ c8 +=1 }
				if($1==9){ c9 +=1 }
				       	print( c0 "\n" c1 "\n" c2 "\n" c3 "\n" c4 "\n" c5 "\n" c6 "\n" c7 "\n" c8 "\n" c9 )}  }' >digit8

cat $1 | awk 'BEGIN{c0=0; c1=0;c2=0;c3=0; c4=0; c5=0; c6=0; c7=0; c8=0; c9=0}
		{ if($2==9){    if($1==0){ c0 +=1 }
				if($1==1){ c1 +=1 }
				if($1==2){ c2 +=1 }
				if($1==3){ c3 +=1 }
				if($1==4){ c4 +=1 }
				if($1==5){ c5 +=1 }
				if($1==6){ c6 +=1 }
				if($1==7){ c7 +=1 }
				if($1==8){ c8 +=1 }
				if($1==9){ c9 +=1 }
				       	print( c0 "\n" c1 "\n" c2 "\n" c3 "\n" c4 "\n" c5 "\n" c6 "\n" c7 "\n" c8 "\n" c9 )}  }' >digit9



for d in $(seq 0 1 9)
	do tail -n 10 digit$d >digitx$d
       	done

paste digitx* 
rm digit*
				