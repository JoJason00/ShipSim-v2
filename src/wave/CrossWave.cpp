#include "WaveBase.h"
#include "CrossWave.h"


CrossWave::CrossWave(const CrossWaveConfig& crosswave)
:config(crosswave){}


double CrossWave::Eta(double t) const
{
	return 0.0;
}
void CrossWave::Exciting(double tn, FKphi& fkphi)
{
}
void CrossWave::loadData(const fkpData& Data)
{
}
double CrossWave::getAmp()
{
	// Representative amplitude for the crossing sea = wave1's amplitude
	// (mirrors getFreq()/direction(), which already delegate to wave1).
	// Feeds the coupled writer's kA / ω_e reference; returning 0 here made
	// every non-dimensional column degenerate for cross seas.
	return config.wave1 ? config.wave1->getAmp() : 0.0;
}
double CrossWave::getFreq()
{
	return config.wave1 ? config.wave1->getFreq() : 0.0;
}
double CrossWave::direction()
{
	return config.wave1 ? config.wave1->direction() : 0.0;
}
