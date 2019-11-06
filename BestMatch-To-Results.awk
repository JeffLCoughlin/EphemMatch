function abs(value){return (value<0?-value:value);}

BEGIN{print("#        TCE       KIC           Source       KIC  Prat    Dist ΔRow ΔCol   Mag2   ΔMag     Drat      Mechanism Flag  Psig  Tsig")}

{
gamma = 0.358153282785417
alpha = 6.93385733756478
sigma = 4.0
stdev = 0.790185331621022


if(NR%3==0)
  {
  type1=$7; name1=$1; kic1=$2; p1=$3; mag1=$8; depth1=$6; pr1=$9; sep=$14; mod1a=$16; out1a=$17; row1a=$18; col1a=$19; mod1b=$21; out1b=$22; row1b=$23; col1b=$24; mod1c=$26; out1c=$27; row1c=$28; col1c=$29; mod1d=$31; out1d=$32; row1d=$33; col1d=$34; psig=$11; tsig=$12;
  };
if(NR%3==1)
  {
  type2=$7; name2=$1; kic2=$2; p2=$3; depth2=$6; mag2=$8; pr2=$9; mod2a=$16; out2a=$17; row2a=$18; col2a=$19; mod2b=$21; out2b=$22; row2b=$23; col2b=$24; mod2c=$26; out2c=$27; row2c=$28; col2c=$29; mod2d=$31; out2d=$32; row2d=$33; col2d=$34;
  };
if(NR%3==1 && NR>2)
  {
  cause="Unknown";
  if(sep<55.0*sqrt(10.0^(-0.4*mag2+6.0)+1.0) || sep<55.0*sqrt(10.0^(-0.4*mag1+6.0)+1.0))
    cause="Direct-PRF";
  if(sep>55.0*sqrt(10.0^(-0.4*mag2+6.0)+1.0) && sep>55.0*sqrt(10.0^(-0.4*mag1+6.0)+1.0) && mod1a!=mod2a && mod1a>0 && mod2a>0)
    cause="Reflection";
  if(sep>55.0*sqrt(10.0^(-0.4*mag2+6.0)+1.0) && sep>55.0*sqrt(10.0^(-0.4*mag1+6.0)+1.0) && ( (abs(row1a-row2a)<=5 && mod1a==mod2a && mod1a>0 && mod2a>0) || (abs(row1b-row2b)<=5 && mod1b==mod2b && mod1b>0 && mod2b>0) || (abs(row1c-row2c)<=5 && mod1c==mod2c && mod1c>0 && mod2c>0) || (abs(row1d-row2d)<=5 && mod1d==mod2d && mod1d>0 && mod2d>0) ) )
    cause="Row-Anomaly";
  if(sep>55.0*sqrt(10.0^(-0.4*mag2+6.0)+1.0) && sep>55.0*sqrt(10.0^(-0.4*mag1+6.0)+1.0) && ( (abs(col1a-col2a)<=5 && mod1a==mod2a && mod1a>0 && mod2a>0) || (abs(col1b-col2b)<=5 && mod1b==mod2b && mod1b>0 && mod2b>0) || (abs(col1c-col2c)<=5 && mod1c==mod2c && mod1c>0 && mod2c>0) || (abs(col1d-col2d)<=5 && mod1d==mod2d && mod1d>0 && mod2d>0) ) )
    cause="Col-Anomaly";
  if(sep>55.0*sqrt(10.0^(-0.4*mag2+6.0)+1.0) && sep>55.0*sqrt(10.0^(-0.4*mag1+6.0)+1.0) && ( (mod1a==mod2a && out1a!=out2a && mod1a>0 && mod2a>0 && out1a>0 && out2a>0 && abs(col1a-col2a)<=5 && abs(row1a-row2a)<=5) || (mod1b==mod2b && out1b!=out2b && mod1b>0 && mod2b>0 && out1b>0 && out2b>0 && abs(col1b-col2b)<=5 && abs(row1b-row2b)<=5) || (mod1c==mod2c && out1c!=out2c && mod1c>0 && mod2c>0 && out1c>0 && out2c>0 && abs(col1c-col2c)<=5 && abs(row1c-row2c)<=5) || (mod1d==mod2d && out1d!=out2d && mod1d>0 && mod2d>0 && out1d>0 && out2d>0 && abs(col1d-col2d)<=5 && abs(row1d-row2d)<=5) ) )
    cause="Cross-Talk";
  if(pr1<1) {pr1=1.0/pr1; pr2=1.0};
  if(pr2<1) {pr2=1.0/pr2; pr1=1.0};
  flag=0;
  flag2=0;
  printflag=0;
  if((cause=="Col-Anomaly" || cause=="Row-Anomaly" || cause=="Reflection") && depth1/depth2 > 0.01 && depth1/depth2 < 100.0)
    flag=1;
  if(cause=="Direct-PRF" && log(depth2/depth1)/log(10.0) < log((10.0^(-0.4*(mag1-mag2)))/(0.5*(gamma^2)/(sep^2+gamma^2) + 0.5*exp(-1.0*(sep/alpha)^2.0)) + 1.0)/log(10.0) - sigma*stdev)
    flag=1;
  if(cause=="Direct-PRF" && log(depth1/depth2)/log(10.0) < log((10.0^(-0.4*(mag2-mag1)))/(0.5*(gamma^2)/(sep^2+gamma^2) + 0.5*exp(-1.0*(sep/alpha)^2.0)) + 1.0)/log(10.0) - sigma*stdev)
    flag2=1;
  if((cause=="Row-Anomaly" || cause=="Col-Anomaly") && out1a!=out2a && out1b!=out2b && out1c!=out2c && out1d!=out2d && flag==0)
    flag=2;
  if((cause=="Row-Anomaly" || cause=="Col-Anomaly") && out1a!=out2a && out1b!=out2b && out1c!=out2c && out1d!=out2d && flag==1)
    flag=3;
  if(depth1<=depth2 || flag2==1)
    {
    printf("%7s %09i %16s %09i %2i:%-2i %7.1f %4i %4i %6.2f %6.2f %9.4E %12s %4i %5.3f %5.3f\n",name1,kic1,name2,kic2,pr2+0.5,pr1+0.5,sep,row1a-row2a,col1a-col2a,mag2,mag1-mag2,depth2/depth1,cause,flag,psig,tsig);          # Added 0.5 to pr1 and pr2 to combat awk's lack of rounding with %i it always rounds down to the nearest integer.
    printflag=1;
    }
  if(flag!=0 && type2=="TCE")
    {
    printf("%7s %09i %16s %09i %2i:%-2i %7.1f %4i %4i %6.2f %6.2f %9.4E %12s %4i %5.3f %5.3f\n",name2,kic2,name1,kic1,pr1+0.5,pr2+0.5,sep,row2a-row1a,col2a-col1a,mag1,mag2-mag1,depth1/depth2,cause,flag,psig,tsig)  # Added 0.5 to pr1 and pr2 to combat awk's lack of rounding with %i it always rounds down to the nearest integer.
    printflag=1;
    }
#  if(printflag==0)
#    printf("FUCK %7s %09i %16s %09i %2i:%-2i %7.1f %4i %4i %6.2f %6.2f %9.4E %12s %4i %5.3f %5.3f\n",name2,kic2,name1,kic1,pr1+0.5,pr2+0.5,sep,row2a-row1a,col2a-col1a,mag1,mag2-mag1,depth1/depth2,cause,flag,psig,tsig)
  }
}

# Flag of 1 = Parent is not likely true parent
# Flag of 2 = Col-Anomaly across different outputs - convoluted but physical
# Flag of 3 = Col-Anom across different outputs and not likely true parent
