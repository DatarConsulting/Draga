#include "gdn-predictive.hpp"
#include "gdn-local.hpp"
#include <array>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <filesystem>
#include <random>
#include <numeric>
#include <sched.h>
namespace fs=std::filesystem;
namespace gp=gdn_predictive;
using Mat=std::array<float,16384>;
using Vec=std::array<float,128>;
using Clock=std::chrono::steady_clock;
constexpr size_t N=16384,HDR=16;
void require(bool v,const char*s){if(!v)throw std::runtime_error(s);}
uint32_t get32(const uint8_t*p){return uint32_t(p[0])|uint32_t(p[1])<<8|uint32_t(p[2])<<16|uint32_t(p[3])<<24;}
void put32(uint8_t*p,uint32_t x){for(int j=0;j<4;j++)p[j]=uint8_t(x>>(8*j));}
struct Spec{std::string name;int kind,rank;std::vector<float> basis;};
// kinds: 0 row magnitude; 1 shared basis; 2 cross approximation;
// 3 Codec A; 4 adaptive Codec B; 5 local; 6 raw.
struct Packet{std::vector<uint8_t>b;size_t model=0,residual=0;unsigned rank=0;double relative_l2=0;};
std::vector<float> load_floats(const fs::path&p){
    std::ifstream f(p,std::ios::binary);require(bool(f),"missing float file");
    auto bytes=fs::file_size(p);require(bytes%4==0,"float size");std::vector<float>v(bytes/4);f.read((char*)v.data(),bytes);require(bool(f),"float read");return v;
}
void prediction(const Spec&s,const uint8_t*meta,size_t nb,unsigned rank,float*p){
    std::fill(p,p+N,0.f);const bool simd=gp::simd_available();
    if(s.kind==0){require(nb==512&&rank==0,"row metadata");
        for(int i=0;i<128;i++){float f;std::memcpy(&f,meta+i*4,4);for(int j=0;j<128;j++)p[i*128+j]=f;}return;}
    require(rank<=96,"rank bound");
    require(nb==rank*(s.kind==1?258:516),"factor metadata");
    if(s.kind==1)require(s.basis.size()==size_t(s.rank)*128&&rank==unsigned(s.rank),"basis dimensions");
    Vec u,v;
    for(unsigned k=0;k<rank;k++){
        if(s.kind==1){std::copy_n(s.basis.data()+k*128,128,u.data());gp::decode_vector(meta+k*258,v.data());}
        else{gp::decode_vector(meta+k*516,u.data());gp::decode_vector(meta+k*516+258,v.data());}
        gp::update(p,nullptr,u.data(),v.data(),simd);
    }
}
Packet encode(const Spec&s,const float*x){
    Packet p;
    if(s.kind>=3){
        p.b.resize(HDR+131072);size_t n=0;
        if(s.kind==3)n=gdn_lossless::encode(x,p.b.data()+HDR);
        if(s.kind==4)n=gp::encode(x,p.b.data()+HDR);
        if(s.kind==5)n=gdn_local::encode(x,p.b.data()+HDR);
        if(s.kind==6){n=65536;std::memcpy(p.b.data()+HDR,x,n);}
        p.b.resize(HDR+n);p.b[0]=0xa7;p.b[1]=s.kind;put32(p.b.data()+12,p.b.size());p.residual=n;return p;
    }
    require(gp::valid_environment(),"floating environment");
    if(s.kind!=0)for(size_t i=0;i<N;i++)require(std::isfinite(x[i]),"nonfinite factor input");
    Mat pred{},remainder;Vec u{},v{};std::vector<uint8_t>meta;const bool simd=gp::simd_available();
    if(s.kind==0){meta.resize(512);
        for(int i=0;i<128;i++){
            // Reference magnitude is a rounded geometric-center word; sign is coded separately.
            uint64_t sum=0;for(int j=0;j<128;j++)sum+=gp::word(x[i*128+j])&0x7fffffff;
            uint32_t ref=uint32_t(sum/128);std::memcpy(meta.data()+4*i,&ref,4);
        }
    }else if(s.kind==1){
        meta.resize(s.rank*258);p.rank=s.rank;
        for(int a=0;a<s.rank;a++){
            std::fill(v.begin(),v.end(),0.f);
            // Fixed FMA order, equal for energy and Hankel bases.
            for(int i=0;i<128;i++)for(int j=0;j<128;j++)v[j]=std::fma(s.basis[a*128+i],x[i*128+j],v[j]);
            for(float value:v)require(std::isfinite(value),"nonfinite basis coefficient");
            require(gp::quantize_vector(v.data(),meta.data()+258*a),"basis coefficient range");
        }
    }else{
        require(s.rank<=32,"cross rank");std::memcpy(remainder.data(),x,65536);meta.resize(s.rank*516);
        for(int a=0;a<s.rank;a++){
            auto index=gp::pivot(remainder.data(),simd);float pivot=remainder[index];if(pivot==0)break;
            auto row=index/128,col=index%128;
            for(int i=0;i<128;i++){u[i]=remainder[i*128+col]/pivot;v[i]=remainder[row*128+i];}
            for(int i=0;i<128;i++)require(std::isfinite(u[i])&&std::isfinite(v[i]),"nonfinite cross coefficient");
            require(gp::quantize_vector(u.data(),meta.data()+a*516)&&gp::quantize_vector(v.data(),meta.data()+a*516+258),"cross coefficients");
            gp::update(pred.data(),remainder.data(),u.data(),v.data(),simd);++p.rank;
        }
        meta.resize(p.rank*516);
    }
    prediction(s,meta.data(),meta.size(),p.rank,pred.data());
    double err=0,energy=0;std::array<uint32_t,N>z;size_t histogram[33]={};
    for(size_t i=0;i<N;i++){
        double e=double(x[i])-pred[i];err+=e*e;energy+=double(x[i])*x[i];
        int64_t d=int64_t(gp::word(x[i])&0x7fffffff)-int64_t(gp::word(pred[i])&0x7fffffff);
        z[i]=uint32_t(d>=0?2*d:-2*d-1);++histogram[z[i]?32-__builtin_clz(z[i]):0];
    }
    p.relative_l2=energy?std::sqrt(err/energy):0.;
    size_t exceptions=N-histogram[0],best=SIZE_MAX,ne=0;unsigned b=32;
    for(unsigned w=0;w<=32;w++){
        size_t len=gp::residual_bytes(w,exceptions);
        if(len<best){best=len;b=w;ne=exceptions;}
        if(w<32)exceptions-=histogram[w+1];
    }
    p.model=meta.size();p.residual=best;p.b.resize(HDR+meta.size()+best,0);
    auto out=p.b.data();out[0]=0xa7;out[1]=s.kind;out[2]=p.rank;out[3]=b;out[4]=ne!=0;gp::put16(out+6,ne);put32(out+8,meta.size());put32(out+12,p.b.size());
    std::copy(meta.begin(),meta.end(),p.b.begin()+HDR);
    auto signs=out+HDR+meta.size(),low=signs+N/8,bitmap=low+N*b/8;
    gp::bit_writer lows{low},highs{bitmap+(ne?N/8:0)};
    uint32_t mask=b==32?0xffffffffu:(uint32_t(1)<<b)-1;
    for(size_t i=0;i<N;i++){
        signs[i/8]|=uint8_t((gp::word(x[i])>>31)<<(i%8));lows.put(z[i]&mask,b);
        if(b<32&&(z[i]>>b)){bitmap[i/8]|=uint8_t(1u<<(i%8));highs.put(z[i]>>b,32-b);}
    }
    lows.finish();highs.finish();require(size_t(highs.p-out)==p.b.size(),"encoded length");return p;
}
void decode(const Spec&s,const Packet&p,float*x){
    const auto in=p.b.data();size_t bytes=p.b.size();require(bytes>=HDR&&in[0]==0xa7&&in[1]==s.kind&&get32(in+12)==bytes,"packet header");
    if(s.kind>=3){
        if(s.kind==3)gdn_lossless::decode(in+HDR,bytes-HDR,x);
        if(s.kind==4)gp::decode(in+HDR,bytes-HDR,x);
        if(s.kind==5)gdn_local::decode(in+HDR,bytes-HDR,x);
        if(s.kind==6){require(bytes==HDR+65536,"raw size");std::memcpy(x,in+HDR,65536);}return;
    }
    unsigned rank=in[2],b=in[3];size_t ne=gp::get16(in+6),nb=get32(in+8);
    require(b<=32&&ne<=N&&in[4]==(ne!=0)&&!in[5]&&!(b==32&&ne),"residual fields");
    require(bytes==HDR+nb+gp::residual_bytes(b,ne),"residual size");
    prediction(s,in+HDR,nb,rank,x);
    const auto signs=in+HDR+nb,low=signs+N/8,bitmap=low+N*b/8;
    if(ne){size_t count=0;for(size_t i=0;i<N/8;i++)count+=__builtin_popcount(unsigned(bitmap[i]));require(count==ne,"exception count");}
    gp::bit_reader lows{low},highs{bitmap+(ne?N/8:0)};
    for(size_t i=0;i<N;i++){
        uint32_t z=lows.get(b);if(ne&&((bitmap[i/8]>>(i%8))&1))z|=highs.get(32-b)<<b;
        int64_t d=(z&1)?-int64_t(z/2)-1:int64_t(z/2);
        int64_t a=int64_t(gp::word(x[i])&0x7fffffff)+d;require(a>=0&&a<=0x7fffffff,"magnitude range");
        uint32_t w=uint32_t(a)|(uint32_t((signs[i/8]>>(i%8))&1)<<31);std::memcpy(x+i,&w,4);
    }
}
void step(float*s,const float*q,const float*k,const float*v,float g,float beta,float*y){
    float decay=std::exp(g),scale=1.f/std::sqrt(128.f);
    for(int j=0;j<128;j++){
        float dot=0;for(int i=0;i<128;i++){s[j*128+i]*=decay;dot+=s[j*128+i]*k[i];}
        float correction=(v[j]-dot)*beta,out=0;
        for(int i=0;i<128;i++){s[j*128+i]+=k[i]*correction;out+=s[j*128+i]*(q[i]*scale);}y[j]=out;
    }
}
volatile uint64_t sink=0;
template<class F>double timer(F f){auto t=Clock::now();f();return std::chrono::duration<double,std::micro>(Clock::now()-t).count();}
struct Record{std::string doc;int offset,head;Mat x;Vec q,k,v;float g,beta;bool update;};
std::vector<Record> records(const fs::path&root,std::ofstream&replay){
    std::vector<Record>out;
    for(auto doc:std::vector<std::string>{"atelier","coast","ledger"}){
        auto base=root/"inputs"/doc;auto state=load_floats(base/"native_state.f32");
        auto q=load_floats(base/"q.f32"),k=load_floats(base/"k.f32"),v=load_floats(base/"v.f32"),g=load_floats(base/"g.f32"),b=load_floats(base/"beta.f32"),observed=load_floats(base/"readout.f32");
        require(state.size()==16*N&&q.size()==33*16*128&&k.size()==q.size()&&v.size()==q.size()&&g.size()==33*16&&b.size()==g.size(),"input shape");
        for(int t=0;t<=32;t++){
            if(t%8==0)for(int h=0;h<16;h++){
                Record r;r.doc=doc;r.offset=t;r.head=h;std::copy_n(state.data()+h*N,N,r.x.data());r.update=t<32;
                int next=std::min(t+1,32);std::copy_n(q.data()+(next*16+h)*128,128,r.q.data());std::copy_n(k.data()+(next*16+h)*128,128,r.k.data());std::copy_n(v.data()+(next*16+h)*128,128,r.v.data());r.g=g[next*16+h];r.beta=b[next*16+h];out.push_back(r);
            }
            if(t==32)break;
            double err=0,den=0,maxerr=0;int next=t+1;
            for(int h=0;h<16;h++){
                Vec y;step(state.data()+h*N,q.data()+(next*16+h)*128,k.data()+(next*16+h)*128,v.data()+(next*16+h)*128,g[next*16+h],b[next*16+h],y.data());
                for(int j=0;j<128;j++){double ref=observed[(next*16+h)*128+j],e=y[j]-ref;err+=e*e;den+=ref*ref;maxerr=std::max(maxerr,std::fabs(e));}
            }
            replay<<doc<<','<<next<<','<<std::sqrt(err/den)<<','<<maxerr<<'\n';
        }
    }return out;
}
void save_frame(const fs::path&path,const Spec&s,const std::vector<Packet>&ps,const std::vector<Mat>&original){
    require(ps.size()==16,"frame heads");std::vector<uint8_t>frame(128,0);std::memcpy(frame.data(),"HANKEL01",8);put32(frame.data()+8,s.kind);put32(frame.data()+12,s.rank);put32(frame.data()+16,16);put32(frame.data()+20,s.basis.size()*4);
    for(int h=0;h<16;h++)put32(frame.data()+32+4*h,ps[h].b.size());
    const auto bp=reinterpret_cast<const uint8_t*>(s.basis.data());if(!s.basis.empty())frame.insert(frame.end(),bp,bp+s.basis.size()*4);
    for(auto&p:ps)frame.insert(frame.end(),p.b.begin(),p.b.end());
    std::ofstream f(path,std::ios::binary);f.write((char*)frame.data(),frame.size());f.close();
    // Reload a self-contained serialized frame; reconstruct basis from its bytes.
    std::ifstream in(path,std::ios::binary);std::vector<uint8_t>back(frame.size());in.read((char*)back.data(),back.size());require(bool(in),"frame read");
    require(!std::memcmp(back.data(),"HANKEL01",8),"frame magic");Spec deployed=s;deployed.basis.resize(get32(back.data()+20)/4);
    if(!deployed.basis.empty())std::memcpy(deployed.basis.data(),back.data()+128,deployed.basis.size()*4);
    size_t pos=128+deployed.basis.size()*4;
    for(int h=0;h<16;h++){auto len=get32(back.data()+32+4*h);require(pos+len<=back.size(),"frame bounds");Packet p;p.b.assign(back.begin()+pos,back.begin()+pos+len);Mat dec;decode(deployed,p,dec.data());require(!std::memcmp(dec.data(),original[h].data(),65536),"saved frame roundtrip");pos+=len;}
    require(pos==back.size(),"frame trailing");
}
void unit_tests(){
    std::mt19937 rng(161106);std::vector<Spec> specs={{"row",0,0,{}},{"identity",1,1,std::vector<float>(128,0.f)},{"cross",2,1,{}}};specs[1].basis[0]=1;
    size_t checks=0;
    for(auto&s:specs)for(int test=0;test<12;test++){
        Mat x{};for(size_t i=0;i<N;i++)x[i]=test==0?0.f:test==1?-(0.f):float(int(rng()%100000)-50000)/1048576.f;
        auto p=encode(s,x.data());Mat y;decode(s,p,y.data());require(!std::memcmp(x.data(),y.data(),65536),"unit roundtrip");++checks;
        auto bad=p;bad.b.pop_back();bool caught=false;try{decode(s,bad,y.data());}catch(...){caught=true;}require(caught,"truncation rejection");++checks;
        bad=p;bad.b[3]=33;caught=false;try{decode(s,bad,y.data());}catch(...){caught=true;}require(caught,"width rejection");++checks;
    }
    // All IEEE words tested using row predictor, including NaN payloads and signed zero.
    Spec row{"row",0,0,{}};for(int trial=0;trial<8;trial++){
        Mat x;for(size_t i=0;i<N;i++){uint32_t w=rng();std::memcpy(x.data()+i,&w,4);}auto p=encode(row,x.data());Mat y;decode(row,p,y.data());require(!std::memcmp(x.data(),y.data(),65536),"random words");++checks;
    }
    std::cout<<"{\"unit_checks\":"<<checks<<",\"passed\":true}\n";
}
int main(int argc,char**argv){
    try{
        require(argc>=2,"experiment ROOT [REPEATS|--test]");if(argc>2&&std::string(argv[2])=="--test"){unit_tests();return 0;}
        fs::path root=argv[1];int repeats=argc>2?std::stoi(argv[2]):3;
        cpu_set_t available;CPU_ZERO(&available);require(!sched_getaffinity(0,sizeof(available),&available),"get affinity");int cpu=0;while(!CPU_ISSET(cpu,&available))++cpu;cpu_set_t chosen;CPU_ZERO(&chosen);CPU_SET(cpu,&chosen);require(!sched_setaffinity(0,sizeof(chosen),&chosen),"pin CPU");
        fs::create_directories(root/"results");fs::create_directories(root/"packets");
        std::ofstream replay(root/"results/replay.csv");replay<<std::setprecision(12)<<"doc,offset,relative_readout_error,max_absolute_error\n";auto data=records(root,replay);
        std::vector<Spec>specs={{"row_ref",0,0,{}},{"codec_A",3,0,{}},{"codec_B",4,0,{}},{"local",5,0,{}},{"raw",6,0,{}}};
        for(int r:{1,2,4,8,16,32,64,96}){
            for(auto label:{"energy","hankel"})specs.push_back({std::string(label)+"_"+std::to_string(r),1,r,load_floats(root/"bases"/(std::string(label)+"_"+std::to_string(r)+".f32"))});
            if(r<=32)specs.push_back({"cross_"+std::to_string(r),2,r,{}});
        }
        std::ofstream sizes(root/"results/sizes.csv"),times(root/"results/times.csv");
        sizes<<std::setprecision(12)<<"doc,offset,head,method,rank,actual_rank,packet_bytes,model_bytes,residual_bytes,basis_bytes,relative_l2,codec_A_bytes\n";
        times<<std::setprecision(12)<<"doc,offset,head,method,trial,encode_us,decode_us,cycle_us\n";
        size_t roundtrips=0,transitions=0,frames=0;std::vector<std::vector<Packet>> saved(specs.size());std::vector<Mat>saved_original;
        for(size_t ri=0;ri<data.size();ri++){
            auto&r=data[ri];auto a=encode(specs[1],r.x.data());
            std::vector<Packet>packets(specs.size());Mat ref=r.x;Vec refy;
            if(r.update)step(ref.data(),r.q.data(),r.k.data(),r.v.data(),r.g,r.beta,refy.data());
            for(size_t j=0;j<specs.size();j++){
                auto&s=specs[j];auto p=encode(s,r.x.data());Mat decoded;decode(s,p,decoded.data());require(!std::memcmp(decoded.data(),r.x.data(),65536),"state roundtrip");++roundtrips;
                if(r.update){Vec y;step(decoded.data(),r.q.data(),r.k.data(),r.v.data(),r.g,r.beta,y.data());require(!std::memcmp(decoded.data(),ref.data(),65536)&&!std::memcmp(y.data(),refy.data(),512),"transition equality");auto next=encode(s,decoded.data());decode(s,next,decoded.data());require(!std::memcmp(decoded.data(),ref.data(),65536),"updated roundtrip");++transitions;}
                sizes<<r.doc<<','<<r.offset<<','<<r.head<<','<<s.name<<','<<s.rank<<','<<p.rank<<','<<p.b.size()<<','<<p.model<<','<<p.residual<<','<<s.basis.size()*4<<','<<p.relative_l2<<','<<a.b.size()<<'\n';
                if(ri<16&&(s.name=="energy_96"||s.name=="hankel_96"||s.name=="cross_32"||s.name=="codec_B"||s.name=="codec_A"||s.name=="row_ref"))saved[j].push_back(p);
                packets[j]=std::move(p);
            }
            if(ri<16)saved_original.push_back(r.x);
            for(int trial=-1;trial<repeats;trial++)for(size_t order=0;order<specs.size();order++){
                size_t j=(order+ri+size_t(trial+1)*7)%specs.size();auto&s=specs[j];Mat x;
                double et=timer([&]{auto p=encode(s,r.x.data());sink=sink+p.b.size();});
                double dt=timer([&]{decode(s,packets[j],x.data());sink=sink+gp::word(x[0]);});
                double ct=-1;
                if(r.update)ct=timer([&]{decode(s,packets[j],x.data());Vec y;step(x.data(),r.q.data(),r.k.data(),r.v.data(),r.g,r.beta,y.data());auto p=encode(s,x.data());sink=sink+p.b.size()+gp::word(y[0]);});
                if(trial>=0)times<<r.doc<<','<<r.offset<<','<<r.head<<','<<s.name<<','<<trial<<','<<et<<','<<dt<<','<<ct<<'\n';
            }
            if((ri+1)%16==0){sizes.flush();times.flush();std::cerr<<"completed "<<ri+1<<"/"<<data.size()<<" heads\n";}
        }
        for(size_t j=0;j<specs.size();j++)if(saved[j].size()==16){save_frame(root/"packets"/(specs[j].name+"_atelier_native.frame"),specs[j],saved[j],saved_original);++frames;}
        std::ofstream checks(root/"results/checks.json");checks<<"{\"roundtrips\":"<<roundtrips<<",\"exact_transitions\":"<<transitions<<",\"self_contained_layer_frames\":"<<frames<<",\"records\":"<<data.size()<<",\"cpu\":"<<cpu<<",\"simd\":"<<(gp::simd_available()?"true":"false")<<",\"passed\":true}\n";
        std::cout<<"Completed; sink="<<sink<<'\n';
    }catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 2;}
    return 0;
}
