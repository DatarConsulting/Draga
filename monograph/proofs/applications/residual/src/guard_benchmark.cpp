#define main original_experiment_main
#include "experiment.cpp"
#undef main
int main(int argc,char**argv){
 try{
    require(argc==2,"guard ROOT");fs::path root=argv[1];
    cpu_set_t available;CPU_ZERO(&available);require(!sched_getaffinity(0,sizeof(available),&available),"affinity");int cpu=0;while(!CPU_ISSET(cpu,&available))++cpu;cpu_set_t one;CPU_ZERO(&one);CPU_SET(cpu,&one);require(!sched_setaffinity(0,sizeof(one),&one),"pin");
    std::ofstream replay(root/"results/guard_replay.csv");auto data=records(root,replay);
    std::ofstream times(root/"results/guard_times.csv");times<<std::setprecision(12)<<"doc,offset,head,rank,trial,encode_us,decode_us,cycle_us,stored_bytes,selected_hankel\n";
    Spec a{"codec_A",3,0,{}};size_t chosen_h=0,checked=0;
    for(int rank:{2,96}){
      Spec h{"hankel",1,rank,load_floats(root/"bases"/("hankel_"+std::to_string(rank)+".f32"))};
      auto enc=[&](const float*x){auto pa=encode(a,x),ph=encode(h,x);return ph.b.size()<pa.b.size()?ph:pa;};
      auto dec=[&](const Packet&p,float*x){decode(p.b[1]==1?h:a,p,x);};
      for(auto&r:data){
        auto p=enc(r.x.data());chosen_h+=p.b[1]==1;Mat x;dec(p,x.data());require(!std::memcmp(x.data(),r.x.data(),65536),"guard roundtrip");++checked;
        if(r.update){Mat ref=r.x;Vec yr,y;step(ref.data(),r.q.data(),r.k.data(),r.v.data(),r.g,r.beta,yr.data());step(x.data(),r.q.data(),r.k.data(),r.v.data(),r.g,r.beta,y.data());auto next=enc(x.data());dec(next,x.data());require(!std::memcmp(x.data(),ref.data(),65536)&&!std::memcmp(y.data(),yr.data(),512),"guard transition");}
        for(int trial=-1;trial<3;trial++){
          auto et=timer([&]{auto packet=enc(r.x.data());sink=sink+packet.b.size();});
          auto dt=timer([&]{dec(p,x.data());sink=sink+gp::word(x[0]);});double ct=-1;
          if(r.update)ct=timer([&]{dec(p,x.data());Vec y;step(x.data(),r.q.data(),r.k.data(),r.v.data(),r.g,r.beta,y.data());auto packet=enc(x.data());sink=sink+packet.b.size()+gp::word(y[0]);});
          if(trial>=0)times<<r.doc<<','<<r.offset<<','<<r.head<<','<<rank<<','<<trial<<','<<et<<','<<dt<<','<<ct<<','<<p.b.size()<<','<<(p.b[1]==1)<<'\n';
        }
      }
    }
    std::ofstream j(root/"results/guard_checks.json");j<<"{\"roundtrips\":"<<checked<<",\"selected_hankel\":"<<chosen_h<<",\"passed\":true,\"protocol\":\"separate sequential guard followup; ranks2 and96; all decoding dispatches by serialized kind\"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
}
