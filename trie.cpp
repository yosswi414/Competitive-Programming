/*
*   Date    : Sep 10 2019
*   Author  : yosswi414
*   Note    : If you find any bugs or some poor points pls tell me, I really want to know that. Thank you.
*/

#include<bits/stdc++.h>
using namespace std;
#define NIL -1

struct node{
  char c;
  int cld,sib,cnt;
  node(char ch,int a=-1,int b=-1):c(ch),cld(a),sib(b),cnt(0){}
};

struct trie{
  vector<node> dat;
  int ptr=0;
  public:
  trie(int siz):dat(vector<node>(siz,node(0))),ptr(0){}
  void addstr(string);
  void addchild(int,int);
  bool contains(string);
  string mss();
  int msc();
};

/*
*   trie::addchild(int par,int child) - make child as one of the children of par
*   It's just a subroutine for trie::addchild() below this one
*/
void trie::addchild(int par,int child){
  vector<node>& dat=this->dat;
  if(dat[par].cld==NIL)dat[par].cld=child;
  else{
    int cur_T=dat[par].cld;
    while(dat[cur_T].sib!=NIL)cur_T=dat[cur_T].sib;
    dat[cur_T].sib=child;
  }
}

/*
*   trie::addstr(string str) - make str belong to trie
*   This part was a big deal(relatively).
*/
void trie::addstr(string str){
  int par=0,cur_T,cur_str=0;
  dat[0].cnt++;
  while(cur_str<(int)str.size()){
    cur_T=dat[par].cld;
    while(cur_T!=NIL&&dat[cur_T].sib!=NIL&&dat[cur_T].c!=str[cur_str]){
      cur_T=dat[cur_T].sib;
    }
    if(cur_T==NIL||dat[cur_T].c!=str[cur_str]){
      while(cur_str<(int)str.size()){
        this->addchild(par,++ptr);
        dat[ptr].c=str[cur_str++];
        dat[ptr].cnt++;
        par=ptr;
      }
    }
    else{
      dat[cur_T].cnt++;
      par=cur_T;
      cur_T=dat[par].cld;
      cur_str++;
    }
  }
}

/*
*   trie::contains(string str) - whether str is contained or not
*   That's it.
*/
bool trie::contains(string str){
  int par=0,cur_T=dat[par].cld,cur_str=0;
  if(str.size()==0)return this->msc();
  while(cur_str<(int)str.size()){
    while(cur_T!=NIL&&dat[cur_T].sib!=NIL&&dat[cur_T].c!=str[cur_str]){
      cur_T=dat[cur_T].sib;
    }
    if(cur_T==NIL||dat[cur_T].c!=str[cur_str])return false;
    else{
      par=cur_T;
      cur_T=dat[par].cld;
      cur_str++;
    }
  }
  return true;
}

int trie::msc(){
  return dat[0].cnt;
}

/*
*   trie::mss() - Most Significant String
*   I just meant "common part" or "intersectional" string, I know I have to change its name with better one
*/

string trie::mss(){
  string res;
  int msc=dat[0].cnt,par=0,cur_T=dat[par].cld;
  while(true){
    while(cur_T!=NIL&&dat[cur_T].sib!=NIL&&dat[cur_T].cnt<msc)cur_T=dat[cur_T].sib;
    if(cur_T==NIL||dat[cur_T].cnt<msc)return res;
    else{
      res+=dat[cur_T].c;
      par=cur_T;
      cur_T=dat[par].cld;
    }
  }
}

/*
*   Usage example problem URL (Japanese only): https://atcoder.jp/contests/yahoo-procon2017-qual/tasks/yahoo_procon2017_qual_c
*/
int main(){
  int n,k;
  cin>>n>>k;
  vector<int> a(k);
  vector<string> s(n);
  for(auto &i:a)cin>>i,i--;
  for(auto &i:s)cin>>i;
  sort(a.begin(),a.end());
  
  trie sub(100000),cmp(100000);
  for(int i=0,t=0;i<n;++i){
    if(t<k&&i==a[t])t++,sub.addstr(s[i]);
    else cmp.addstr(s[i]);
  }
  
  if(!cmp.contains(sub.mss())){
    int t=sub.mss().size();
    string res=sub.mss();
    while(t>0&&!cmp.contains(res.substr(0,t-1)))t--;
    cout<<res.substr(0,t);
  }
  else cout<<-1;
  return 0;
}
