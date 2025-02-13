# Reads Hao's orbital sim output
# c0 cov - day - descend
# c1 cov - day - ascend
# c2 cov - night - descend
# c3 cov - night - ascend
# c4 pow - day   - descend
# c5 pow - day   - ascend
# c6 pow - night - descend
# c7 pow - night - ascend

BEGIN{
  nCases=8;  # number of cases on the input
  srand(int(seed));
  numb=0;
}


($1!="id"){  # read data
  for(i=0;i<nCases;i++){
    nIn[numb,i]=$(i+4);
  }
  numb++;
}

END{
  # turn data into a PDF of ascend and descend for now
  maxIn=0;
  for(i=0;i<numb;i++){
    n[i,0]=n[i,1]=0;
    for(j=0;j<nCases;j++){
      k=j%2;
      n[i,k]+=nIn[i,j];
    }
    if(n[i,0]>maxIn)maxIn=n[i,0];
    if(n[i,1]>maxIn)maxIn=n[i,1];
  }

  # bin up histogram
  for(i=0;i<=maxIn;i++){
    for(j=0;j<=maxIn;j++){
      hist[i,j]=0;
    }
  }
  total=0;
  for(i=0;i<numb;i++){
    # apply cloud fraction

    nA=nD=0;
    for(j=0;j<n[i,0];j++)if(rand()>=cFrac)nA++;
    for(j=0;j<n[i,1];j++)if(rand()>=cFrac)nD++;

    hist[nA,nD]++;
    total++;
  }

  # write out PDF
  for(i=0;i<=maxIn;i++){
    for(j=0;j<=maxIn;j++){
      hist[i,j]/=total;
      print i,j,hist[i,j];
    }
  }
}

