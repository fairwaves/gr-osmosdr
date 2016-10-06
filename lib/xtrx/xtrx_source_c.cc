/* -*- c++ -*- */
/*
 * Copyright 2012 Dimitri Stolnikov <horiz0n@gmx.net>
 *
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#include <fstream>
#include <string>
#include <sstream>
#include <map>

#include <boost/assign.hpp>
#include <boost/algorithm/string.hpp>

#include <gnuradio/io_signature.h>
#include <gnuradio/blocks/deinterleave.h>
#include <gnuradio/blocks/float_to_complex.h>

#include "xtrx_source_c.h"

#include "arg_helpers.h"

using namespace boost::assign;


xtrx_source_c_sptr make_xtrx_source_c(const std::string &args)
{
  return gnuradio::get_initial_sptr(new xtrx_source_c(args));
}

static size_t parse_nchan(const std::string &args)
{
  size_t nchan = 1;

  dict_t dict = params_to_dict(args);

  if (dict.count("nchan"))
    nchan = boost::lexical_cast< size_t >( dict["nchan"] );

  if (nchan < 1)
    nchan = 1;
  else if (nchan > 2)
    nchan = 2;

  return nchan;
}

xtrx_source_c::xtrx_source_c(const std::string &args) :
  gr::sync_block("xtrx_source_c",
                 gr::io_signature::make(0, 0, 0),
                 gr::io_signature::make(parse_nchan(args),
                                        parse_nchan(args),
                                        sizeof(gr_complex))),
  _xtrxdev(NULL),
  _rate(0),
  _master(0),
  _freq(0),
  _corr(0),
  _bandwidth(0),
  _auto_gain(false),
  _otw(XTRX_WF_16),
  _mimo_mode(parse_nchan(args) > 1),
  _rema(false)
{

  dict_t dict = params_to_dict(args);

  if (dict.count("otw_format")) {
    const std::string& otw = dict["otw_format"];
    if (otw == "sc16" || otw == "16") {
      _otw = XTRX_WF_16;
    } else if (otw == "sc12" || otw == "12") {
      _otw = XTRX_WF_12;
    } else if (otw == "sc8" || otw == "8") {
      _otw = XTRX_WF_8;
    } else {
      throw std::runtime_error("Parameter `otw_format` should be {sc16,sc12,sc8}");
    }
  }

  if (dict.count("master")) {
    _master = boost::lexical_cast< double >( dict["master"]);
  }

  _channels = parse_nchan(args);

/*
  if (dict.count("direct_samp"))
    direct_samp = boost::lexical_cast< unsigned int >( dict["direct_samp"] );

  if (dict.count("offset_tune"))
    offset_tune = boost::lexical_cast< unsigned int >( dict["offset_tune"] );
*/

  std::cerr << args.c_str() << std::endl;

  int loglevel = 4;
  if (dict.count("loglevel")) {
    loglevel = boost::lexical_cast< int >( dict["loglevel"] );
  }

  bool lmsreset = 0;
  if (dict.count("lmsreset")) {
    lmsreset = boost::lexical_cast< bool >( dict["lmsreset"] );
  }

  unsigned xtrxflag = (loglevel & XTRX_O_LOGLVL_MASK) | ((lmsreset) ? XTRX_O_RESET : 0);
  std::cerr << "xtrx_source_c::xtrxflag = " << xtrxflag << std::endl;
  int res = xtrx_open("/dev/xtrx0", xtrxflag, &_xtrxdev);
  if (res) {
    std::stringstream message;
    message << "Couldn't open "  ": Error: " << -res;
    throw std::runtime_error( message.str() );
  }

  std::cerr << "xtrx_source_c::xtrx_source_c()" << std::endl;
}

xtrx_source_c::~xtrx_source_c()
{
  std::cerr << "xtrx_source_c::~xtrx_source_c()" << std::endl;

  if (_xtrxdev)
    xtrx_close(_xtrxdev);
}

std::string xtrx_source_c::name()
{
  return "GrLibXTRX";
}

std::vector<std::string> xtrx_source_c::get_devices( bool fake )
{
  std::vector<std::string> devices;

  // TODO
  devices += "/dev/xtrx0";

  return devices;
}

size_t xtrx_source_c::get_num_channels( void )
{
  return output_signature()->max_streams();
}

osmosdr::meta_range_t xtrx_source_c::get_sample_rates( void )
{
  osmosdr::meta_range_t range;
  range += osmosdr::range_t( 1000000, 160000000, 1 );
  return range;
}

double xtrx_source_c::set_sample_rate( double rate )
{
  std::cerr << "Set sample rate " << rate << std::endl;

  int res = xtrx_set_samplerate(_xtrxdev, _master, rate, 0,
                                NULL, &_rate, NULL);
  if (res) {
    std::cerr << "Unable to set samplerate, error=" << res << std::endl;
  }
  return get_sample_rate();
}

double xtrx_source_c::get_sample_rate( void )
{
  return _rate;
}

osmosdr::freq_range_t xtrx_source_c::get_freq_range( size_t chan )
{
  osmosdr::freq_range_t range;
  range += osmosdr::range_t( double(0.1e9), double(3.8e9), 1); // as far as we know
  return range;
}

double xtrx_source_c::set_center_freq( double freq, size_t chan )
{
  _freq = freq;
  double corr_freq = (freq)*(1.0 + (_corr) * 0.000001);

  std::cerr << "Set freq " << freq << std::endl;

  int res = xtrx_tune(_xtrxdev, XTRX_TUNE_RX_FDD, corr_freq, &_freq);
  if (res) {
    std::cerr << "Unable to deliver frequency " << corr_freq << std::endl;
  }

  return get_center_freq(chan);
}

double xtrx_source_c::get_center_freq( size_t chan )
{
  return _freq;
}

double xtrx_source_c::set_freq_corr( double ppm, size_t chan )
{
  _corr = ppm;

  set_center_freq(_freq, chan);

  return get_freq_corr( chan );
}

double xtrx_source_c::get_freq_corr( size_t chan )
{
  return _corr;
}

static const std::map<std::string, xtrx_gain_type_t> s_lna_map = boost::assign::map_list_of
    ("LNA", XTRX_RX_LNA_GAIN)
    ("TIA", XTRX_RX_TIA_GAIN)
    ("PGA", XTRX_RX_PGA_GAIN)
;

static xtrx_gain_type_t get_gain_type(const std::string& name)
{
  std::map<std::string, xtrx_gain_type_t>::const_iterator it;

  it = s_lna_map.find(name);
  if (it != s_lna_map.end()) {
    return it->second;
  }

  return XTRX_RX_LNA_GAIN;
}

static const std::vector<std::string> s_lna_list = boost::assign::list_of
    ("LNA")("TIA")("PGA")
;

std::vector<std::string> xtrx_source_c::get_gain_names( size_t chan )
{
  return s_lna_list;
}

osmosdr::gain_range_t xtrx_source_c::get_gain_range( size_t chan )
{
  return get_gain_range("LNA", chan);
}

osmosdr::gain_range_t xtrx_source_c::get_gain_range( const std::string & name, size_t chan )
{
  osmosdr::gain_range_t range;

  if (name == "LNA") {
    range += osmosdr::range_t( 0, 24,  3 );
    range += osmosdr::range_t( 25, 30, 1 );
  } else if (name == "TIA") {
    range += osmosdr::range_t( 0 );
    range += osmosdr::range_t( 9 );
    range += osmosdr::range_t( 12 );
  } else if (name == "PGA") {
    range += osmosdr::range_t( -12.5, 12.5, 1 );
  }

  return range;
}

bool xtrx_source_c::set_gain_mode( bool automatic, size_t chan )
{
  _auto_gain = automatic;
  return get_gain_mode(chan);
}

bool xtrx_source_c::get_gain_mode( size_t chan )
{
  return _auto_gain;
}

double xtrx_source_c::set_gain( double gain, size_t chan )
{
  return set_gain(gain, "LNA", chan);
}

double xtrx_source_c::set_gain( double igain, const std::string & name, size_t chan )
{
  osmosdr::gain_range_t gains = xtrx_source_c::get_gain_range( name, chan );
  double gain = gains.clip(igain);
  double actual_gain;
  xtrx_gain_type_t gt = get_gain_type(name);

  std::cerr << "Set gain " << name << " (" << gt << "): " << igain << std::endl;

  int res = xtrx_set_gain(_xtrxdev, /*(chan == 0) ? XTRX_CH_A : XTRX_CH_B*/ XTRX_CH_AB, gt, gain, &actual_gain);
  if (res) {
    std::cerr << "Unable to set gain `" << name.c_str() << "`; err=" << res << std::endl;
  }

  switch (gt) {
  case XTRX_RX_LNA_GAIN: _gain_lna = actual_gain; break;
  case XTRX_RX_TIA_GAIN: _gain_tia = actual_gain; break;
  case XTRX_RX_PGA_GAIN: _gain_pga = actual_gain; break;
  default: break;
  }

  return actual_gain;
}

double xtrx_source_c::get_gain( size_t chan )
{
  return get_gain("LNA");
}

double xtrx_source_c::get_gain( const std::string & name, size_t chan )
{
  xtrx_gain_type_t gt = get_gain_type(name);
  switch (gt) {
  case XTRX_RX_LNA_GAIN: return _gain_lna;
  case XTRX_RX_TIA_GAIN: return _gain_tia;
  case XTRX_RX_PGA_GAIN: return _gain_pga;
  }
  return 0;
}

double xtrx_source_c::set_if_gain(double gain, size_t chan)
{
  return set_gain(gain, "PGA", chan);
}

double xtrx_source_c::set_bandwidth( double bandwidth, size_t chan )
{
  std::cerr << "Set bandwidth " << bandwidth << " chan " << chan << std::endl;

  if (bandwidth <= 0.0) {
      bandwidth = get_sample_rate() * 0.75;
      if (bandwidth < 1.4e6) {
          bandwidth = 1.4e6;
      }
  }

  int res = xtrx_tune_rx_bandwidth(_xtrxdev, (chan == 0) ? XTRX_CH_A : XTRX_CH_B, bandwidth, &_bandwidth);
  if (res) {
    std::cerr << "Can't set bandwidth: " << res << std::endl;
  }
  return get_bandwidth(chan);
}

double xtrx_source_c::get_bandwidth( size_t chan )
{
  return _bandwidth;
}


static const std::map<std::string, xtrx_antenna_t> s_ant_map = boost::assign::map_list_of
    ("RXW", XTRX_RX_W)
    ("RXL", XTRX_RX_L)
    ("RXH", XTRX_RX_H)
;
static const std::map<xtrx_antenna_t, std::string> s_ant_map_r = boost::assign::map_list_of
    (XTRX_RX_W, "RXW")
    (XTRX_RX_L, "RXL")
    (XTRX_RX_H, "RXH")
;

static xtrx_antenna_t get_ant_type(const std::string& name)
{
  std::map<std::string, xtrx_antenna_t>::const_iterator it;

  it = s_ant_map.find(name);
  if (it != s_ant_map.end()) {
    return it->second;
  }

  return XTRX_RX_W;
}

static const std::vector<std::string> s_ant_list = boost::assign::list_of
    ("RXW")("RXL")("RXH")
;


std::vector< std::string > xtrx_source_c::get_antennas( size_t chan )
{
  return s_ant_list;
}

std::string xtrx_source_c::set_antenna( const std::string & antenna, size_t chan )
{
  _ant = get_ant_type(antenna);

  std::cerr << "Set antenna " << antenna << std::endl;

  int res = xtrx_set_antenna(_xtrxdev, _ant);
  if (res) {
    std::cerr << "Can't set antenna: " << antenna << std::endl;
  }
  return get_antenna( chan );
}

std::string xtrx_source_c::get_antenna( size_t chan )
{
  return s_ant_map_r.find(_ant)->second;
}

int xtrx_source_c::work (int noutput_items,
                         gr_vector_const_void_star &input_items,
                         gr_vector_void_star &output_items)
{
  //std::cerr << "work(" << noutput_items << ")" << "\n";
  static float dummy[65536];

  float* outdata = (float*)output_items[0];
  float* outdata2 = dummy;

  if (output_items.size() > 1) {
      outdata2 =  (float*)output_items[1];
  }

  int res = xtrx_recv_sync(_xtrxdev, noutput_items * 2, outdata, outdata2);
  if (res) {
    std::stringstream message;
    message << "xtrx_recv_sync error: " << -res;
    throw std::runtime_error( message.str() );
  }

  return noutput_items;
}

bool xtrx_source_c::start()
{
  //TODO:
  std::cerr << "xtrx_source_c::start(otw=" << _otw << ")" << std::endl;
  int res = xtrx_run(_xtrxdev, XTRX_RX, _otw, /*(_channels == 1) ? XTRX_CH_A :*/ XTRX_CH_AB , XTRX_IQ_FLOAT32);
  if (res) {
    std::cerr << "Got error: " << res << std::endl;
  }

  _rema = false;
  return res == 0;
}

bool xtrx_source_c::stop()
{
  //TODO:
  std::cerr << "xtrx_source_c::stop()" << std::endl;
  int res = xtrx_stop(_xtrxdev, XTRX_RX);
  if (res) {
    std::cerr << "Got error: " << res << std::endl;
  }

  return res == 0;
}
