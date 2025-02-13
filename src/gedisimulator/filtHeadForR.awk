BEGIN{
  first=1;
}

($0){
  if($1=="#"){
    for(i=2;i<=NF;i++){
      notNumb=1;
      for(j=0;j<1000;j++){
        if($i==j){
          notNumb=0;
          break;
        }
      }
      if(notNumb==1){
        printf(" %s",$i);
      }
    }
    printf("\n");
  }
}

