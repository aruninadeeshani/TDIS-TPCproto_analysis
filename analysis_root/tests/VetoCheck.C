#include "../TPCDetector.h"
#include <chrono>
// Targeted fixture: no raw inputs, geometry maps or plot refresh needed.
void VetoCheck() {
    using namespace PersonalTPC;
    auto folder=fs::temp_directory_path()/("tpc-veto-check-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(folder/"source/auxiliary");
    fs::create_directories(folder/"results/voltage");
    struct Cleanup {fs::path path;~Cleanup(){fs::remove_all(path);}} cleanup{folder};
    auto write=[](const fs::path& path,const std::string& text){std::ofstream f(path);f<<text;};
    auto db=folder/"cuts.ini";
    write(db,"[global]\nsetup=large_medium\n[run 1605]\ndrift_voltage_v=\n");
    Config config(folder.c_str(),(folder/"results").c_str(),"",db.c_str());
    config.source=folder/"source";
    auto bundled=config.Auxiliary("rejected_events.csv");
    write(bundled,"run,event\n1605,1054\n1605,1055\n1616,7\n");
    auto live=config.output/"voltage/rejected_events.csv",summary=live.parent_path()/"run_summary.csv";
    Check(config.VetoFile(1605)==bundled,"Absent live ledger did not fall back");
    write(live,"run,event\n1605,9999\n1606,8\n");
    Check(config.VetoFile(1605)==bundled,"Ledger without processing summary superseded snapshot");
    write(summary,"run,noise_events,selected\n1606,1,False\n1616,0,True\n");
    Check(config.VetoFile(1605)==bundled && ReadVetoes(1605,config)==std::set<int>({1054,1055}),
          "Unknown-voltage skipped run lost its two bundled vetoes");
    Check(config.VetoFile(1616)==live && ReadVetoes(1616,config).empty(),
          "Processed clean run inherited stale bundled vetoes");
    Check(config.VetoFile(1606)==live && ReadVetoes(1606,config)==std::set<int>({8}),
          "Processed but unselected run ignored live vetoes");
    auto explicitLedger=folder/"explicit.csv";write(explicitLedger,"run,event\n1605,42\n");
    config.vetoOverride=explicitLedger;
    Check(config.VetoFile(1605)==explicitLedger && ReadVetoes(1605,config)==std::set<int>({42}),
          "Explicit veto override lost precedence");
    config.vetoOverride.clear();
    write(live,"run,event\n");write(summary,"noise_events,run\n0,1605\n");
    Check(config.VetoFile(1605)==live && ReadVetoes(1605,config).empty(),
          "Processed run with empty ledger incorrectly fell back");
    write(summary,"unrelated,column\n0,1605\n");
    Check(config.VetoFile(1605)==bundled,"Summary without run column authorized live ledger");
    auto mixed=config.output/"voltage/large_medium/rejected_events.csv";
    fs::create_directories(mixed.parent_path());
    write(mixed,"run,event\n1605,99\n");write(mixed.parent_path()/"run_summary.csv","run\n1616\n");
    Check(config.VetoFile(1605)==bundled,"Mixed ledger bypassed run processing check");
    write(mixed.parent_path()/"run_summary.csv","run\n1605\n");
    Check(config.VetoFile(1605)==mixed && ReadVetoes(1605,config)==std::set<int>({99}),
          "Mixed ledger ignored its processed run");
    Config ms=config;ms.outputOverride=false;ms.setupOverride="medium_small";ms.runInfo.clear();
    auto msLive=config.output/"setups/medium_small/voltage/rejected_events.csv";
    fs::create_directories(msLive.parent_path());write(msLive,"run,event\n");
    write(msLive.parent_path()/"run_summary.csv","run\n1616\n");
    Check(ms.VetoFile(1616)==msLive && ReadVetoes(1616,ms).empty() && ms.VetoFile(1605)==bundled,
          "MS ledger processing check depends on a missing map or leaks across runs");
    std::cout<<"PASS: native veto coverage gate, preserved skipped-run vetoes, processed zero-veto runs, explicit overrides and mixed/MS paths."<<std::endl;
}
