#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "FWCore/Utilities/interface/Exception.h"

const HcalTrigTowerDetId HcalTrigTowerDetId::Undefined(0x4a000000u);

HcalTrigTowerDetId::HcalTrigTowerDetId() {
}


HcalTrigTowerDetId::HcalTrigTowerDetId(uint32_t rawid) : DetId(rawid) {
}

HcalTrigTowerDetId::HcalTrigTowerDetId(int ieta, int iphi) : DetId(Hcal,HcalTriggerTower) {
  int depth = 1;  // forcing depth = 1 for consistency with emap;
  int version = 0;
  id_|= ((depth&0x7)<<14) |
        ((ieta>0)?(0x2000|(ieta<<7)):((-ieta)<<7)) |
    (iphi&0x7F);
  id_|=((version&0x7)<<17);
}

HcalTrigTowerDetId::HcalTrigTowerDetId(int ieta, int iphi, int depth) : DetId(Hcal,HcalTriggerTower) {
  const int ones = depth % 10;
  const int tens = (depth - ones) / 10;
  int mdepth = 1, mversion = 0; // default value
  if (tens == 0) {
     mdepth = ones;
     mversion = 0;
  }else if(tens == 1){
     mdepth = 0;
     mversion = 1;
  }else{
     assert(tens <= 1);
  }
  id_|=((mdepth&0x7)<<14) |
    ((ieta>0)?(0x2000|(ieta<<7)):((-ieta)<<7)) |
    (iphi&0x7F);
  id_|=((mversion&0x7)<<17);
}

HcalTrigTowerDetId::HcalTrigTowerDetId(int ieta, int iphi, int depth, int version) : DetId(Hcal,HcalTriggerTower) {
  id_|=((depth&0x7)<<14) |
    ((ieta>0)?(0x2000|(ieta<<7)):((-ieta)<<7)) |
    (iphi&0x7F);
  id_|=((version&0x7)<<17);
}
 
HcalTrigTowerDetId::HcalTrigTowerDetId(const DetId& gen) {
  if (!gen.null() && (gen.det()!=Hcal || gen.subdetId()!=HcalTriggerTower)) {
    throw cms::Exception("Invalid DetId") << "Cannot initialize HcalTrigTowerDetId from " << std::hex << gen.rawId() << std::dec; 
  }
  id_=gen.rawId();
}

void HcalTrigTowerDetId::setVersion(int version) {
  id_|=((version&0x7)<<17);
}

HcalTrigTowerDetId& HcalTrigTowerDetId::operator=(const DetId& gen) {
  if (!gen.null() && (gen.det()!=Hcal || gen.subdetId()!=HcalTriggerTower)) {
    throw cms::Exception("Invalid DetId") << "Cannot assign HcalTrigTowerDetId from " << std::hex << gen.rawId() << std::dec; 
  }
  id_=gen.rawId();
  return *this;
}

std::ostream& operator<<(std::ostream& s,const HcalTrigTowerDetId& id) {
  s << "(HcalTrigTower v" << id.version() << ": " << id.ieta() << ',' << id.iphi();
  if (id.depth()>0) s << ',' << id.depth();
  
  return s << ')';
}


