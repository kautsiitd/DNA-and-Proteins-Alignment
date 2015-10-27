#include<bits/stdc++.h>

using namespace std;

vector<vector<int> > MC;
int CC;
map<vector<int>, int> MCost;
int size_V,K;
vector<char> V;
vector<string> XStr;
vector<string> CurrentStrings;
map<string, int> CostMap;
vector<vector<int> > XInd, XRemain;
int NODES_expanded;
vector<int> Lengths;
map<pair<int, int>, vector<vector<int> > > HCostMap;
vector<vector<vector<vector<int> > > > HMatrices;

bitset<28> BitStart;
vector<string> XSoFar;
int AnswerDepth,MinNextCost;

bool IDAstarDeeperBit(int Costbound, int depth,int CurrCost,vector<int> CurrIndices,vector<int> RemLengths)
;

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


vector<int> ConvertStringToVec(string Str,map<char, int> MyMap)
//convert the 'n' length string to an 'n' length vector
{
    vector<int> Ind(Str.length());
    for(int i=0;i<Str.length();i++)
    {
        Ind[i]=MyMap[Str[i]];           //converts string to vec<int>
    }
    return Ind;
}

int MinCostDashes(vector<int> RemainLen)
//finds the minimum number of dashes required to form equal length strings
{
    
    int maxlen = RemainLen[0] , minDashes=0;
    
    for(int i=1;i<K;i++)
    {
        if(RemainLen[i]>maxlen)
        {
            maxlen=RemainLen[i];                        //compute maxlen
        }
    }
    for(int i=0;i<K;i++)
    {
        minDashes+= maxlen-RemainLen[i];                //get min. dashes
    }
    
    return minDashes*CC;
}


void GenerateHeuristicCost()
//one time execution. Creates Pairwise Cost for all the possible substrings and stores in matrix, refered to by indices of the pair, ie for strings(vectors) i,j, the cost matrix is HMAtrices[i-1][j]
{
    vector<int> v1,v2;
    int t1,t2,t3,t4,t5,t6;
    int L1,L2;
    pair<int, int> MyPair;
    for(int i=1;i<K;i++)
    {
        vector<vector<vector<int> > > IMatrices;
        for(int j=0;j<i;j++)
        {
            v1=XInd[i];
            v2=XInd[j];
            reverse(v1.begin(), v1.end());              //reverse strings for ease of future reference
            reverse(v2.begin(), v2.end());
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
            IMatrices.push_back(Matrix);
            
        }
        HMatrices.push_back(IMatrices);
    }
}


int HeuristicCost(vector<int> RemLen)
//uses matrix cost and min. dash cost to find heuristic cost
{
    int HCost=0;
    pair<int, int> MyPair;
    
    vector<vector<int > > Mat;
    for(int i=1;i<K;i++)
    {
        for(int j=0;j<i;j++)
        {
            HCost += HMatrices[i-1][j][RemLen[i]][RemLen[j]];       //add PairWise Costs
        }
    }
    HCost += MinCostDashes(RemLen);                                 //add costs for Dashes
    return HCost;
    
}


bool IDAstarBit(int Costbound, int depth, bitset<28> bitnum,int CurrCost,vector<int> CurrIndices,vector<int> RemLengths)
//Given a bitset number, makes the child and calculate Cost. if Cost<Bound, goes deeper,else returns false

{
    vector<int> child(K);
    string childStr;
    
    //generate the vector representing this child
    for(int i=0;i<K;i++)
    {
        if(bitnum.test(i))                          //next Child has char
        {
            //add dash if allcharacters ended
            if(RemLengths[i]!=0)
            {
                RemLengths[i]--;
                childStr+=XStr[i][CurrIndices[i]];
                
                CurrIndices[i]++;
            }
            else
            {
                childStr+='-';
            }
        }
        else
        {
            childStr+='-';
        }
    }
    //end of child generation
    
    //calc total(incl. hueristic) cost for that child

    int OSC= CostMap[childStr];
    
    int ActualCost= OSC+ CurrCost;
    
    int HCOST = HeuristicCost(RemLengths);                //Change Parameter
    
    int Cost= ActualCost + HCOST ;              //lower bound for this cost
    
    if(Cost<=Costbound)
    {
        XSoFar[depth]=childStr;
        //XCurr.push_back(child);
        return IDAstarDeeperBit(Costbound, depth+1, ActualCost, CurrIndices, RemLengths);    //child is transpose!!!!!
    }
    else
    {
        //if not found, generate next child
        if (Cost<MinNextCost)
        {
            MinNextCost= Cost;             //do only if cost>costbound??. ensure its the deepest level.
        }
    }
    return false;
}



bool IDAstarDeeperBit(int Costbound, int depth,int CurrCost,vector<int> CurrIndices,vector<int> RemLengths)
//generates all nodes at one depth. Exits if answer is found
{
   
    NODES_expanded++;
    //generate one child
    bitset<28> bitnum= BitStart;
    bool Found=false;
    bool AllEmpty=true;
    
    //Check if ALL empty
    for(int i=0;i<K;i++)
    {
        if(RemLengths[i]==0)
        {
            bitnum.set(1,1);
        }
        else
        {
            AllEmpty=false;
        }
        
    }

    if(AllEmpty)                    //no string remaining,and no failure=> success
    {
        AnswerDepth = depth;
        return true;
    }
    while(!bitnum.test(K) and !Found)                   //iterate through all possoble nodes at this depth
    {
        NODES_expanded++;
        
        Found = IDAstarBit(Costbound, depth, bitnum, CurrCost, CurrIndices, RemLengths);
        bitnum = AddOneBit(bitnum);
    }
    
    
    return Found;
}


void searchBit()
//called to initiate Search
{
    bool Found=false;
    int bound=1;
    BitStart.reset();
    /*
     for(int i=K;i<28;i++)
     {
     bitStart.flip(i);                       //define the bitnumber,set when states covered
     }
     */
    BitStart.flip(0);                       //bitnum=0 => ALL dashes, not wanted =>increased to 1
    GenerateHeuristicCost();
    int NumIter=0;
    vector<int> InitIndices(K,0);
    vector<vector<int> > empty;
    
    while(!Found)
    {
        NumIter++;
        MinNextCost= numeric_limits<int>::max();
        XRemain=XInd;
        
        
        Found = IDAstarDeeperBit(bound , 0, 0, InitIndices, Lengths);
        //NextBound is the lower bound on actual cost, using min. Estimated Cost in nodes expanded
        cout<<NumIter<<" "<<MinNextCost<<" "<<bound<<"\n";
        
        if(bound>=MinNextCost)                            //as a precaution, to insure increment
        {
            bound++;
        }
        else
        {
            bound=MinNextCost;                //allows faster increment went certain answer cost is higher
        }
    }
    
//    cout<<NumIter<<" Iterations of Cost increase were required \n";
    
}


int main(int argc, const char * argv[])
{
    //cout<<argc<<endl;
    //Take input
    fstream File;
    //string inputName,outputName;
    //cin>>inputName;

    File.open("input.txt",ios::in);
    
    float Time;
    int i=0, j;
    char comm;
    File>>Time>>size_V;
    Time= 60*Time;                      //convert to seconds
    V.resize(size_V);
    map<char, int> MyMap;
    while(i<(size_V-1))
    {
        File>>V[i]>>comm;                //input of vocab
        MyMap[V[i]] = i;                //stores map for each letter
        i++;
    }
    File>>V[i];
    MyMap[V[i]] = i;
    File>>K;
    XStr.resize(K);
    vector<int> temp;
    i=0;
    int totalLengthOfString=0;
    
    while(i<K)
    {
        File>>XStr[i];                   //input of K strings
        temp= ConvertStringToVec(XStr[i],MyMap);
        totalLengthOfString +=XStr[i].length();
        Lengths.push_back(XStr[i].length());
        XInd.push_back(temp);
        i++;
    }
    File>>CC;
    i=0;
    while(i<(size_V+1))
    {
        j=0;
        vector<int> temp(size_V+1);
        while (j<(size_V+1))
        {
            File>>temp[j];
            j++;
        }
        MC.push_back(temp);
        i++;
    }
    char hash;
    File>>hash;
    vector<int> empty;
    
    File.close();
    //End of Input
    
    
    
    NODES_expanded=0;
    CreateCostMap(empty, 0, "");

//    cout<<"Cost Map creation = "<<time_spent<<"\n";

    XSoFar.resize(totalLengthOfString);                     //maximum possible length of Answer
    searchBit();
    File.open("output.txt",ios::out);
    for(int i=0;i<XSoFar[0].size();i++)
    {
        for(int j=0;j<AnswerDepth;j++)
        {
            // string temp = XSoFar[j];
            //cout<<XSoFar[j][i];
            File<<XSoFar[j][i];
        }
        File<<"\n";
        //cout<<"\n";
    }
    File.close();


    
    
    
//    cout<<NODES_expanded<<" nodes expanded\n";
//    cout<<" Hueristic took time = "<<time_spent<<"\nCostCalc took = "<<time_costCalc;
//    cout<<" Presearch = "<<time_presearch<<"\n"<<"Time Calling = "<<time_calling<<"\n";
//    time_spent = (double)(endO - beginO) / CLOCKS_PER_SEC;
//    cout<<"Time Taken ="<<time_spent<<"\n \n";
    
    return 0;
}

