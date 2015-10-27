#include<bits/stdc++.h>
using namespace std;
#define inf 1e9
#define ll long long

// Global Variables
time_t start;

vector<string> CurrentStrings;
map<string, int> CostMap;
vector< int> InitLengths;
vector< vector<int> > HCostMatrix;
map<int, vector<int> > VectorMap;
map<int, vector<vector<int> > > SolnMap;
int HMatSize;
	// input variables
ll nodes=0;
float Time;
int size_V, K;
vector<char> V(28);											// vocab letters
map<char, int> MyMap;										// each letter Index in MC	{ (T 2),.. }
vector<string> XStr(100);									// input strings
int Strlen[100];											// length of each input string
int CC,MC[28][28];											// conversion cost
vector< vector<int> > XInd(100);							// strings converted to integer array
	// process variable
int maxlen=0;												// store maximum number of letter in vector of strings
															// for calculating number of dashes to make equal length strings
vector<string> ans;
ll ans_len;
vector< vector<int> > ansi;
ll upperbound;
int total_letters=0;
int ans_cost_till=0;
vector<string> temps;										// temp for random solution generated
vector< vector<int> > tempi(20);
int templ=0;
	// Random
int Index[100];
int pos_ran[20][100];
	// climb
int dec_cost;
int no_of_swap;
vector<string> final;
int finall;
// variables end here

void CreateCostMap(vector<int> VecTillNow,int Cost, string StrTillNow)


// This calculates the cost of any possible step(as per the Vocab) and stores it in a Map
{
    //increase size of vector, if size=K, put in Map, else Expand further
    int tempCost;
    vector<int> temp;
    for(int i=0;i<size_V;i++)
    {
        
        temp= VecTillNow;
        temp.push_back(i);
        tempCost = Cost;
        for(int j=0; j<VecTillNow.size();j++)
        {
            tempCost= tempCost + MC[i][VecTillNow[j]];          //increase cost for so far
        }
        if(temp.size()<K)                                       //check size
        {
            CreateCostMap(temp, tempCost, StrTillNow+V[i]);                      //increase size
        }
        else
        {
            CostMap[StrTillNow+V[i]] = tempCost;                           //add to map
        }
    }
    //Code for '-'
    temp= VecTillNow;
    temp.push_back(size_V);
    
    tempCost = Cost;
    for(int j=0; j<VecTillNow.size();j++)
    {
        tempCost= tempCost + MC[size_V][VecTillNow[j]];          //increase cost for so far
    }
    tempCost+=CC;
    
    if(temp.size()<K)                                           //check size
    {
        CreateCostMap(temp, tempCost, StrTillNow+'-');
    }
    else
    {
        CostMap[StrTillNow+'-'] = tempCost;
    }
    return;
}


int CostColumn(string Str)
{
    int Cost=0;
    if(Str[0]=='-')
    {
        Cost+=CC;
    }
    for(int i=1;i<Str.length();i++)
    {
        for(int j=0;j<i;j++)
        {
            Cost+= MC[MyMap[Str[i]]][MyMap[Str[j]]];
        }
        if(Str[i]=='-')
        {
            Cost+=CC;
        }
    }
    return Cost;
}



bitset<28> AddOneBit(bitset<28> bitnum)
//increases the bitset representation by one. Seems to run faster like this
{
    bool ToSwitch;
    for(int i=K;i>=0;i--)
    {
        ToSwitch=true;
        int j=0;
        while(ToSwitch && j<i)
        {
            ToSwitch= (ToSwitch && bitnum.test(j));
            j++;
        }
        if(ToSwitch)
        {
            bitnum.flip(i);
        }
    }
    return bitnum;
}

void GenerateHeuristicCost()
//one time execution. Creates Pairwise Cost for all the possible substrings and stores in matrix, refered to by indices of the pair, ie for strings(vectors) i,j, the cost matrix is HMAtrices[i-1][j]
{
    vector<int> v1,v2;
    HCostMatrix.resize(K);
    
    vector<vector<int> > HMatrix(K,vector<int>(K));
    int t1,t2,t3,t4,t5,t6;
    int L1,L2;
    pair<int, int> MyPair;
    for(int i=1;i<K;i++)
    {
        for(int j=0;j<i;j++)
        {
            v1=XInd[i];
            v2=XInd[j];
            L1= v1.size();
            L2=v2.size();
            vector<vector<int> > Matrix(L1+1,vector<int>(L2+1,0));               //check if 0 init.
            for (int k=1; k<=L1; k++)
            {
                Matrix[k][0]= Matrix[k-1][0]+ MC[v1[k-1]][size_V];
            }
            for (int k=1; k<=L2; k++)
            {
                Matrix[0][k]= Matrix[0][k-1]+ MC[size_V][v2[k-1]];
            }
            for(int k=1;k<=L1;k++)
            {
                for(int l=1;l<=L2;l++)
                {
                    t1= Matrix[k-1][l-1];
                    t2= MC[v1[k-1]][v2[l-1]];
                    t3= Matrix[k-1][l] ;
                    t4= MC[v1[k-1]][size_V];
                    t5=  Matrix[k][l-1] ;
                    t6= MC[size_V][v2[l-1]];
                    Matrix[k][l] = min(t1+t2  ,min( t3+t4 , t5+t6));
                }
            }
            HMatrix[i][j] = Matrix[L1][L2];
            HMatrix[j][i] = Matrix[L1][L2];
        }
    }
    HCostMatrix = HMatrix;
    
}

vector<vector<int> > MatchStrings(vector<vector<int> > v1, vector<vector<int> > v2)
{
    vector<vector<int> > Final1 (v1.size());
    vector<vector<int> > Final2 (v2.size());
    int L1= v1[0].size();
    int v1_Size = v1.size();
    int L2=v2[0].size();
    int v2_Size =v2.size();
    for(int m=0;m<v1_Size;m++)
    {
        reverse(v1[m].begin(), v1[m].end());              //reverse strings for ease of future reference
    }
    for(int m=0;m<v2_Size;m++)
    {
        reverse(v2[m].begin(), v2[m].end());              //reverse strings for ease of future reference
    }
    
    int t1,t2,t3,t4,t5,t6;
    vector<vector<int> > Matrix(L1+1,vector<int>(L2+1,0)),  DirMatrix(L1+1,vector<int>(L2+1,1));               //check if 0 init.

    for(int k=1;k<=L1;k++ )
    {
        DirMatrix[k][0] = 2;
    }
    
    for (int k=1; k<=L1; k++)
    {
        Matrix[k][0]= Matrix[k-1][0] + v2_Size*CC;
        for(int m=0;m<v2_Size;m++)
        {
            for(int n=0;n<v1_Size;n++)
            {
                Matrix[k][0]+= MC[v1[n][k-1]][size_V];
            }
        }
        
    }
    for (int k=1; k<=L2; k++)
    {
        Matrix[0][k]= Matrix[0][k-1] + v1_Size*CC;
        for(int m=0;m<v1_Size;m++)
        {
            for(int n=0;n<v2_Size;n++)
            {
                Matrix[0][k]+= MC[size_V][v2[n][k-1]] ;
            }
        }
    }
    
    
    for(int k=1;k<=L1;k++)
    {
        for(int l=1;l<=L2;l++)
        {
            t1= Matrix[k-1][l-1];
            t2=0;
            for(int m=0;m<v1_Size;m++)
            {
                for(int n=0;n<v2_Size;n++)
                {
                    t2+= MC[v1[m][k-1]][v2[n][l-1]];
                }
            }
            t3= Matrix[k-1][l] ;
            t4= v2_Size*CC;
            for(int m=0;m<v2_Size;m++)
            {
                for(int n=0;n<v1_Size;n++)
                {
                    t4+= MC[v1[n][k-1]][size_V];
                }
            }
            
            t5 = Matrix[k][l-1] ;
            t6=v1_Size*CC;
            for(int m=0;m<v1_Size;m++)
            {
                for(int n=0;n<v2_Size;n++)
                {
                    t6+= MC[size_V][v2[n][l-1]] ;
                }
            }
            t1 = t1+t2;
            t2 = t3+t4;
            t3 = t5+t6;
            Matrix[k][l] = min(t1  ,min( t2 , t3));
            
            if(t1<=t2 && t1<=t3)
            {
                DirMatrix[k][l] = 3;
            }
            else if (t2<= t3)
            {
                DirMatrix[k][l] = 2;
            }
        }
    }
    int k=L1, l=L2;
    
    
    while(k>0 && l>0)
    {
        if(DirMatrix[k][l]==3)
        {
            k--;
            l--;
            for(int m=0;m<v1_Size;m++)
            {
                Final1[m].push_back(v1[m][k]);
            }
            for(int m=0;m<v2_Size;m++)
            {
                Final2[m].push_back(v2[m][l]);
            }
           
        }
        else if(DirMatrix[k][l]==2)
        {
            k--;
            for(int m=0;m<v1_Size;m++)
            {
                Final1[m].push_back(v1[m][k]);
            }
            for(int m=0;m<v2_Size;m++)
            {
                Final2[m].push_back(size_V);
            }
        }
        else
        {
            l--;
            for(int m=0;m<v1_Size;m++)
            {
                Final1[m].push_back(size_V);
            }
            for(int m=0;m<v2_Size;m++)
            {
                Final2[m].push_back(v2[m][l]);
            }
        }
    }
    while(k>0)
    {
        k--;
        for(int m=0;m<v1_Size;m++)
        {
            Final1[m].push_back(v1[m][k]);
        }
        for(int m=0;m<v2_Size;m++)
        {
            Final2[m].push_back(size_V);
        }
    }
    while (l>0)
    {
        l--;
        for(int m=0;m<v1_Size;m++)
        {
            Final1[m].push_back(size_V);
        }
        for(int m=0;m<v2_Size;m++)
        {
            Final2[m].push_back(v2[m][l]);
        }
    }
    for(int n=0;n<v2_Size;n++)
    {
        Final1.push_back(Final2[n]);
    }
    return Final1;
}



pair<int, int> MinPair()
{
    pair<int, int> MyPair;
    int min = numeric_limits<int>::max();
    for(int i=1;i<=HMatSize;i++)
    {
        for(int j=0;j<i;j++)
        {
            if(HCostMatrix[i][j]<min)
            {
                min =HCostMatrix[i][j];
                MyPair = make_pair(i, j);
            }
        }
    }
    return MyPair;
}

vector<vector<int> > ClustalW()
{
    GenerateHeuristicCost();
    vector<int> Sum(K);
    int NumStr = K;
    pair<int,int> MyPair;
    
    while(NumStr!=0)
    {
        vector<int> NextVec;
        vector<vector<int> > NewSoln;
        int firstNum, secondNum;
        vector<vector<int> > First, Second;
        
        MyPair = MinPair();
        firstNum = MyPair.second;
        secondNum = MyPair.first;
        //cout<<firstNum<<" "<<secondNum<<" = "<<HMatSize+1<<"\n";
        
        if(firstNum<K)
        {
            First.push_back(XInd[firstNum]);
            NextVec.push_back(firstNum);
        }
        else
        {
            vector<int> temp = VectorMap[firstNum];
            for(int i=0;i<temp.size();i++)
            {
                vector<vector<int> > tempMat = SolnMap[firstNum];
                vector<int> tempVec = tempMat[i];
                First.push_back(tempVec);
                NextVec.push_back(temp[i]);
            }
            
        }
        if(secondNum<K)
        {
            Second.push_back(XInd[secondNum]);
            NextVec.push_back(secondNum);
        }
        else
        {
            vector<int> temp = VectorMap[secondNum];
            for(int i=0;i<temp.size();i++)
            {
                vector<vector<int> > tempMat = SolnMap[secondNum];
                vector<int> tempVec = tempMat[i];
                Second.push_back(tempVec);
                NextVec.push_back(temp[i]);
            }
        }
        NewSoln =  MatchStrings(First, Second);
        HMatSize++;
        VectorMap[HMatSize] = NextVec;
        SolnMap[HMatSize] = NewSoln;
        
        vector<int> BestVec= HCostMatrix[firstNum];
        for(int i=0;i<HMatSize;i++)
        {
            if(HCostMatrix[secondNum][i]<BestVec[i])
            {
                BestVec[i] = HCostMatrix[secondNum][i];
            }
            HCostMatrix[i].push_back(BestVec[i]);
        }
        BestVec.push_back(numeric_limits<int>::max());
        HCostMatrix.push_back(BestVec);
        for(int i=0;i<NextVec.size();i++)
        {
            for(int j=0;j<=HMatSize;j++)
            {
                HCostMatrix[NextVec[i]][j] = numeric_limits<int>::max();
                HCostMatrix[j][NextVec[i]] = numeric_limits<int>::max();
            }
        }
        if(firstNum>=K)
        {
            for(int j=0;j<=HMatSize;j++)
            {
                HCostMatrix[firstNum][j] = numeric_limits<int>::max();
                HCostMatrix[j][firstNum] = numeric_limits<int>::max();
            }
        }
        if(secondNum>=K)
        {
            for(int j=0;j<=HMatSize;j++)
            {
                HCostMatrix[secondNum][j] = numeric_limits<int>::max();
                HCostMatrix[j][secondNum] = numeric_limits<int>::max();
            }
        }
        
        if(NextVec.size()>=K) {
			vector<vector<int> > ReturnAns(K);
			for(int i=0;i<K;i++) {
				ReturnAns[NextVec[i]] = NewSoln[i];
			}
			return ReturnAns;
		}
        NumStr = K-NextVec.size();
    }
    return XInd;    //NEVER CALLED
}

// Convert string to integer labelling
vector<int> ConvertStringToVec(string Str){
    vector<int> Ind(Str.length());
    for(int i=0;i<Str.length();i++){
        Ind[i]=MyMap[Str[i]];
    }
    return Ind;
}

// finding total_cost (type=0 => for ans type=1 => for temp)
	// checked
ll total_cost(int type,int len){
	ll cost=0;
	for(int i=0;i<len;i++){
		for(int j=0;j<K;j++){
			for(int k=j+1;k<K;k++){
				if(type==0){
					//cout<<ansi[j][i]<<" "<<ansi[k][i]<<endl;
					cost+=MC[ansi[j][i]][ansi[k][i]];
				}
				else if(type==2){
					cost+=MC[MyMap[final[j][i]]][MyMap[final[k][i]]];
				}
				else{
					cost+=MC[tempi[j][i]][tempi[k][i]];
				}
			}
		}
		if(type==0){
			for(int j=0;j<K;j++){
				if(ansi[j][i]==size_V){
					cost+=CC;
				}
			}
		}
		else if(type==2){
			for(int j=0;j<K;j++){
				if(final[j][i]=='-'){
					cost+=CC;
				}
			}
		}
		else{
			for(int j=0;j<K;j++){
				if(tempi[j][i]==size_V){
					cost+=CC;
				}
			}
		}
	}
	return cost;
}

// resetting all Index values to zero
void reset_Index(){
	for(int i=0;i<K;i++){
		Index[i]=0;
	}
}

// generating random solution and storing in temps and tempi
	// checked
void Random(){
	srand (time(NULL));
	// initially Index is 0 in XStr[i]
	reset_Index();
	temps.clear();
	tempi.clear();
	templ=0;
	// temps contains empty strings initially
	for(int i=0;i<K;i++){
		temps.push_back("");
		
	}
	ll letter_remain=total_letters;
	while(letter_remain>0){
		//cout<<endl;
		vector<int> temp;
		string tempstr="";
		int should=0;												// to exclude only - columns
		for(int i=0;i<K;i++){
			int Prob=rand()%100+0;
			int rnd=Prob%2;
			//cout<<Index[i]<<" ";
			if(Index[i]==Strlen[i] || rnd==0){
				tempstr+='-';
				temp.push_back(size_V);
			}
			else{
				tempstr+=XStr[i][Index[i]];
				temp.push_back(XInd[i][Index[i]]);
				letter_remain-=1;
				should=1;
			}
		}
		//cout<<letter_remain<<endl;
		if(should==1){
			//cout<<tempstr<<endl;
			for(int i=0;i<K;i++){
				if(tempstr[i]!='-'){
					// updating pos_ran
					pos_ran[i][Index[i]]=templ;
					Index[i]+=1;
				}
				temps[i]+=tempstr[i];
				tempi[i].push_back(temp[i]);
			}
			templ+=1;
		}
	}
}

void climb(){
	no_of_swap=0;
	dec_cost=0;
	for(int i=0;i<templ;i++){
		for(int j=0;j<K;j++){
			if(temps[j][i]=='-' && Index[j]<Strlen[j]){
				//cout<<"go inside "<<temps[j][i]<<" "<<j<<" "<<i<<endl;
				// find change in cost due to swapping
					// current cost
				ll current_column_cost=0;
				ll after_swap_cost=0;
				for(int k=0;k<K;k++){
					// '-' matching with others
					current_column_cost+=MC[tempi[k][i]][size_V];
					// first letter in row, with others in that's column
					// other..... first letter
					current_column_cost+=MC[tempi[k][pos_ran[j][Index[j]]]][tempi[j][pos_ran[j][Index[j]]]];
					// ignoring CC cost because that would not be part of change in cost
					// '-' matching with others in column of first letter
					after_swap_cost+=MC[tempi[k][pos_ran[j][Index[j]]]][size_V];
					// first letter in row, with others in '-'s column
					after_swap_cost+=MC[tempi[k][i]][tempi[j][pos_ran[j][Index[j]]]];
				}
				// should swap if cost is decreasing
				if(current_column_cost>=after_swap_cost-2*MC[tempi[j][pos_ran[j][Index[j]]]][size_V]){
					// swap change in temps and tempi
					// update pos_ran
					// update dec_cost
					// update no_of_swap
					// update Index[i]
					//cout<<current_column_cost<<" "<<after_swap_cost<<endl;
					//cout<<"swaping "<<"in string "<<j<<" and in column "<<i<<" and "<<pos_ran[j][Index[j]]<<" "<<current_column_cost-after_swap_cost<<endl;
					temps[j][i]=XStr[j][Index[j]];
					temps[j][pos_ran[j][Index[j]]]='-';
					tempi[j][i]=XInd[j][Index[j]];
					tempi[j][pos_ran[j][Index[j]]]=size_V;
					pos_ran[j][Index[j]]=i;
					dec_cost+=current_column_cost-after_swap_cost+2*MC[tempi[j][pos_ran[j][Index[j]]]][size_V];
					no_of_swap+=1;
					Index[j]+=1;
					/*for(int i=0;i<K;i++){
				    	cout<<temps[i]<<endl;
				    	/*for(int j=0;j<templ;j++){
				    		cout<<pos_ran[i][j]<<" ";
				    		//cout<<tempi[j][i]<<" ";
				    	}
				    	//cout<<endl;
				    }
				    cout<<endl<<endl;*/
				}
			}
			else{
				Index[j]+=1;
			}
		}
	}
	//cout<<"swap "<<no_of_swap<<endl;
}

void change_ans(int type,int length){
	ans.clear();
	ans_len=length;
	if(type==0){
		for(int i=0;i<K;i++){
			ans.push_back(temps[i]);
		}
	}
	else{
		for(int i=0;i<K;i++){
			ans.push_back(final[i]);
		}
	}
}

int main(int argc, const char * argv[]){
	time_t start=time(NULL);
	fstream File;
    File.open(argv[1],ios::in);
	upperbound=inf;
    // input start
    char comm;
    File>>Time;
    File>>size_V;
    Time=Time*60;
    
    // vocab input and labelling of letter
    for(int i=0;i<size_V;i++){
    	string data;
    	File>>data;
        V[i]=data[0];
        MyMap[V[i]]=i;
    }
    MyMap['-']=size_V;
    char tempCharForNow = V[size_V];
    V[size_V] = '-';
    // string input and conversion to integer according to labelling and storing length of strings
    File>>K;
    for(int i=0;i<K;i++){
        File>>XStr[i];
    	XInd[i]=ConvertStringToVec(XStr[i]);
    	Strlen[i]=XStr[i].length();
    	total_letters+=Strlen[i];
    	if(maxlen<Strlen[i]){
    		maxlen=Strlen[i];
    	}
    }
    
    // conversion cost input
    File>>CC;
    for(int i=0;i<size_V+1;i++){
    	ll min=inf;
        for(int j=0;j<size_V+1;j++){
        	int data;
        	File>>data;
            MC[i][j]=data;
            if(min>data && data!=0){
            	min=data;
            }
        }
        MC[i][size_V+1]=min;
    }
    File.close();
    // end of input
    // Pranit program
    HMatSize = K-1;
    vector<int> empty;
    //CreateCostMap(empty, 0, "");
    //Search(Time);
   
    vector<vector<int> > Answer = ClustalW();
    ans.clear();
    ansi=Answer;
    for(int i=0;i<Answer.size();i++)
    {
    	string lol="";
        for(int j =0;j<Answer[i].size();j++)
        {
            lol+=V[Answer[i][j]];
        }
        ans.push_back(lol);
    }
    upperbound=total_cost(0,ansi[0].size());
    V[size_V] = tempCharForNow;
    // first dash extended solution
    	// checked
    for(int i=0;i<K;i++){
    	vector<int> tempi;
    	string temps=XStr[i];
    	for(int j=0;j<Strlen[i];j++){
    		tempi.push_back(XInd[i][j]);
    	}
    	for(int j=0;j<maxlen-Strlen[i];j++){
    		temps+='-';
    		tempi.push_back(size_V);
    	}
    	ansi.push_back(tempi);
    	ans.push_back(temps);
    }
    ans_len=maxlen;
    // finding total cost of ans as sending type=0 and len=maxlen;
    if(upperbound>total_cost(0,maxlen)){
    	upperbound=total_cost(0,maxlen);
    }
    //cout<<upperbound<<endl;
    // generating first random solution
    time_t current=time(NULL);
    while((double)(current-start)+5<Time){
		//cout<<upperbound<<endl;
	    Random();
	    ll random_cost=total_cost(1,templ);
	    //cout<<random_cost<<endl;
	    if(random_cost<upperbound){
	    	upperbound=random_cost;
	    	change_ans(0,templ);
	    }
	    reset_Index();
	    dec_cost=inf;
	    no_of_swap=inf;
	    time_t mid=time(NULL);
	    while((double)(mid-start)+5<Time && no_of_swap!=0){
	    	climb();
	    	mid=time(NULL);
	    }
	    finall=0;
	    final.clear();
	    for(int i=0;i<K;i++){
	    	final.push_back("");
	    }
	    for(int i=0;i<templ;i++){
	    	int should=0;
	    	for(int j=0;j<K;j++){
	    		if(temps[j][i]!='-'){
	    			should=1;
	    			break;
	    		}
	    	}
	    	if(should==1){
	    		for(int j=0;j<K;j++){
	    			final[j]+=temps[j][i];
	    		}
	    		finall+=1;
	    	}
	    }
	    ll after_cost=total_cost(2,finall);
	    //cout<<after_cost<<endl;
	    if(after_cost<upperbound){
	    	upperbound=after_cost;
	    	change_ans(1,finall);
	    }
	    current=time(NULL);
	}
	File.open(argv[2], ios::out);
	for(int i=0;i<K;i++){
		File<<ans[i]<<endl;
	}
	File.close();
	//cout<<upperbound<<endl;
	//clock_t end=time(NULL);
	//cout<<(double)(end-start)<<endl;
    return 0;
}
