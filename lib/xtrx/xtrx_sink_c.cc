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

#include "xtrx_sink_c.h"

#include "arg_helpers.h"

using namespace boost::assign;

xtrx_sink_c_sptr make_xtrx_sink_c(const std::string &args)
{
  return gnuradio::get_initial_sptr(new xtrx_sink_c(args));
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

xtrx_sink_c::xtrx_sink_c(const std::string &args) :
  gr::sync_block("xtrx_sink_c",
                 gr::io_signature::make(parse_nchan(args),
                                        parse_nchan(args),
                                        sizeof(gr_complex)),
                 gr::io_signature::make(0, 0, 0)),
  _xtrxdev(NULL),
  _rate(0),
  _master(0),
  _freq(0),
  _corr(0),
  _bandwidth(0),
  _auto_gain(false),
  _otw(XTRX_WF_16),
  _mimo_mode(parse_nchan(args) > 1),
  _rema(false),
  _ts(8192)
{

  dict_t dict = params_to_dict(args);
/*  Not supported yet!
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
*/
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
  std::cerr << "xtrx_sink_c::xtrxflag = " << xtrxflag << std::endl;
  int res = xtrx_open("/dev/xtrx0", xtrxflag, &_xtrxdev);
  if (res) {
    std::stringstream message;
    message << "Couldn't open "  ": Error: " << -res;
    throw std::runtime_error( message.str() );
  }

  if (dict.count("refclk")) {
	xtrx_set_ref_clk(_xtrxdev, boost::lexical_cast< unsigned >( dict["refclk"] ), XTRX_CLKSRC_INT);
  }

  std::cerr << "xtrx_sink_c::xtrx_sink_c()" << std::endl;
  set_alignment(32);
  set_output_multiple(1024);
}

xtrx_sink_c::~xtrx_sink_c()
{
  std::cerr << "xtrx_sink_c::~xtrx_sink_c()" << std::endl;
  if (_xtrxdev) {
    xtrx_close(_xtrxdev);
    _xtrxdev = NULL;
  }
}

std::string xtrx_sink_c::name()
{
  return "GrLibXTRX";
}

std::vector<std::string> xtrx_sink_c::get_devices( bool fake )
{
  std::vector<std::string> devices;

  // TODO
  devices += "/dev/xtrx0";

  return devices;
}

size_t xtrx_sink_c::get_num_channels( void )
{
  return input_signature()->max_streams();
}

osmosdr::meta_range_t xtrx_sink_c::get_sample_rates( void )
{
  osmosdr::meta_range_t range;
  range += osmosdr::range_t( 1000000, 160000000, 1 );
  return range;
}

double xtrx_sink_c::set_sample_rate( double rate )
{
  std::cerr << "Set sample rate " << rate << std::endl;

  int res = xtrx_set_samplerate(_xtrxdev, _master, 0, rate,
                                NULL, NULL, &_rate);
  if (res) {
    std::cerr << "Unable to set samplerate, error=" << res << std::endl;
  }
  return get_sample_rate();
}

double xtrx_sink_c::get_sample_rate( void )
{
  return _rate;
}

osmosdr::freq_range_t xtrx_sink_c::get_freq_range( size_t chan )
{
  osmosdr::freq_range_t range;
  range += osmosdr::range_t( double(0.1e9), double(3.8e9), 1); // as far as we know
  return range;
}

double xtrx_sink_c::set_center_freq( double freq, size_t chan )
{
  _freq = freq;
  double corr_freq = (freq)*(1.0 + (_corr) * 0.000001);

  std::cerr << "Set freq " << freq << std::endl;

  int res = xtrx_tune(_xtrxdev, XTRX_TUNE_TX_FDD, corr_freq, &_freq);
  if (res) {
    std::cerr << "Unable to deliver frequency " << corr_freq << std::endl;
  }

  return get_center_freq(chan);
}

double xtrx_sink_c::get_center_freq( size_t chan )
{
  return _freq;
}

double xtrx_sink_c::set_freq_corr( double ppm, size_t chan )
{
  _corr = ppm;

  set_center_freq(_freq, chan);

  return get_freq_corr( chan );
}

double xtrx_sink_c::get_freq_corr( size_t chan )
{
  return _corr;
}


static const std::vector<std::string> s_lna_list = boost::assign::list_of
    ("TX")
;

std::vector<std::string> xtrx_sink_c::get_gain_names( size_t chan )
{
  return s_lna_list;
}

osmosdr::gain_range_t xtrx_sink_c::get_gain_range( size_t chan )
{
  return get_gain_range("TX", chan);
}

osmosdr::gain_range_t xtrx_sink_c::get_gain_range( const std::string & name, size_t chan )
{
  osmosdr::gain_range_t range;
  range += osmosdr::range_t( -10 );
  range += osmosdr::range_t( 7 );
  range += osmosdr::range_t( 25 );
  return range;
}

bool xtrx_sink_c::set_gain_mode( bool automatic, size_t chan )
{
  _auto_gain = automatic;
  return get_gain_mode(chan);
}

bool xtrx_sink_c::get_gain_mode( size_t chan )
{
  return _auto_gain;
}

double xtrx_sink_c::set_gain( double gain, size_t chan )
{
  return set_gain(gain, "TX", chan);
}

double xtrx_sink_c::set_gain( double igain, const std::string & name, size_t chan )
{
  osmosdr::gain_range_t gains = xtrx_sink_c::get_gain_range( name, chan );
  double gain = gains.clip(igain);
  double actual_gain;

  std::cerr << "Set TX gain: " << igain << std::endl;

  int res = xtrx_set_gain(_xtrxdev, /*(chan == 0) ? XTRX_CH_A : XTRX_CH_B*/ XTRX_CH_AB,
                          XTRX_TX_PAD_GAIN, gain, &actual_gain);
  if (res) {
    std::cerr << "Unable to set gain `" << name.c_str() << "`; err=" << res << std::endl;
  }

  _gain_tx = actual_gain;
  return actual_gain;
}

double xtrx_sink_c::get_gain( size_t chan )
{
  return get_gain("TX");
}

double xtrx_sink_c::get_gain( const std::string & name, size_t chan )
{
  return _gain_tx;
}

double xtrx_sink_c::set_bandwidth( double bandwidth, size_t chan )
{
  std::cerr << "Set bandwidth " << bandwidth << " chan " << chan << std::endl;

  if (bandwidth <= 0.0) {
      bandwidth = get_sample_rate() * 0.75;
      if (bandwidth < 0.5e6) {
          bandwidth = 0.5e6;
      }
  }

  int res = xtrx_tune_tx_bandwidth(_xtrxdev, (chan == 0) ? XTRX_CH_A : XTRX_CH_B, bandwidth, &_bandwidth);
  if (res) {
    std::cerr << "Can't set bandwidth: " << res << std::endl;
  }
  return get_bandwidth(chan);
}

double xtrx_sink_c::get_bandwidth( size_t chan )
{
  return _bandwidth;
}


static const std::map<std::string, xtrx_antenna_t> s_ant_map = boost::assign::map_list_of
    ("B1", XTRX_TX_L)
    ("B2", XTRX_TX_W)
;
static const std::map<xtrx_antenna_t, std::string> s_ant_map_r = boost::assign::map_list_of
    (XTRX_TX_L, "B1")
    (XTRX_TX_W, "B2")
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
    ("B1")("B2")
;


std::vector< std::string > xtrx_sink_c::get_antennas( size_t chan )
{
  return s_ant_list;
}

std::string xtrx_sink_c::set_antenna( const std::string & antenna, size_t chan )
{
  _ant = get_ant_type(antenna);

  std::cerr << "Set antenna " << antenna << std::endl;

  int res = xtrx_set_antenna(_xtrxdev, _ant);
  if (res) {
    std::cerr << "Can't set antenna: " << antenna << std::endl;
  }
  return get_antenna( chan );
}

std::string xtrx_sink_c::get_antenna( size_t chan )
{
  return s_ant_map_r.find(_ant)->second;
}

static const float dummy[65536] = {0.0};

int xtrx_sink_c::work (int noutput_items,
                         gr_vector_const_void_star &input_items,
                         gr_vector_void_star &output_items)
{
  //std::cerr << "tx_work(" << noutput_items << ")" << "\n";
  const float* indata = (const float*)input_items[0];
  const float* indata2 = dummy;

  if (input_items.size() > 1) {
      indata2 =  (float*)input_items[1];
  }

  // Hack send using small bursts
  int remaining = noutput_items;
  unsigned off = 0;
  do {
    int burstsz = remaining;
    if (burstsz > 1024) {
      burstsz = 1024;
    }
    remaining -= burstsz;

    //std::cerr << "tx_work(REM=" << remaining << " TS=" << _ts << ")" << std::endl;
    int res = xtrx_send_burst_sync(_xtrxdev, _ts, 2*burstsz, indata2 + 2*off, indata + 2*off);
    if (res) {
      std::cerr << "Err: " << res << std::endl;

      std::stringstream message;
      message << "xtrx_send_burst_sync error: " << -res;
      throw std::runtime_error( message.str() );
    }

    _ts += burstsz;
    off += burstsz;
  } while (remaining > 0);
  //_ts

  consume(0, noutput_items);
  if (input_items.size() > 1) {
    consume(1, noutput_items);
  }
  return 0;
}

bool xtrx_sink_c::start()
{
  //TODO:
  std::cerr << "xtrx_sink_c::start(otw=" << _otw << ")" << std::endl;
  int res = xtrx_run(_xtrxdev, XTRX_TX, _otw, /*(_channels == 1) ? XTRX_CH_A :*/ XTRX_CH_AB , XTRX_IQ_FLOAT32, 0, 0);
  if (res) {
    std::cerr << "Got error: " << res << std::endl;
  }

  _rema = false;
  return res == 0;
}

bool xtrx_sink_c::stop()
{
  //TODO:
  std::cerr << "xtrx_sink_c::stop()" << std::endl;
  int res = xtrx_stop(_xtrxdev, XTRX_TX);
  if (res) {
    std::cerr << "Got error: " << res << std::endl;
  }

  return res == 0;
}
