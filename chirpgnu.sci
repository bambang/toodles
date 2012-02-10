function y = chirp(t, f0, t1, f1, form, phase)
//Evaluate a chirp signal at time t.  A chirp signal is a frequency swept cosine wave.
// Calling Sequence
//  y = chirp(t)
//  y = chirp(t,f0)
//  y = chirp(t,f0,t1)
//  y = chirp(t,f0,t1,f1)
//  y = chirp(t,f0,t1,f1,form)
//  y = chirp(t,f0,t1,f1,form,phase)
// Parameters
// t: vector of times to evaluate the chirp signal
// f0: frequency at time t=0 [0 Hz]
// t1: time t1 [1 sec]
// f1: frequency at time t=t1 [100 Hz]
// form: shape of frequency sweep;   'linear' :     f(t) = (f1-f0)*(t/t1) + f0, ,    'quadratic':   f(t) = (f1-f0)*(t/t1)^2 + f0,    'logarithmic': f(t) = (f1-f0)^(t/t1) + f0
// phase: phase shift at t=0
// Description
// If you want a different sweep shape f(t), use the following:
//    y = cos(2*%pi*integral(f(t)) + 2*%pi*f0*t + phase);
// Examples
//     tfrsp(chirp([0:0.001:5])',1:5001,128,'plot'); // linear, 0-100Hz in 1 sec
//     tfrsp(chirp([-2:0.001:15], 400, 10, 100, 'quadratic')',1:5001,128,'plot');
//    plot(chirp([0:1/8000:5], 200, 2, 500, "logarithmic"));
//
// // Shows linear sweep of 100 Hz/sec starting at zero for 5 sec
// // since the sample rate is 1000 Hz, this should be a diagonal
// // from bottom left to top right.
//  tfrsp(chirp([0:0.001:5])',1:5001,128,'plot'); // linear, 0-100Hz in 1 sec
//
// // Shows a quadratic chirp of 400 Hz at t=0 and 100 Hz at t=10
// // Time goes from -2 to 15 seconds.
// stacksize('max');
// tfrsp(chirp([-2:0.001:15], 400, 10, 100, 'quadratic')',1:17001,128,'plot');
//
// // Shows a logarithmic chirp of 200 Hz at t=0 and 500 Hz at t=2
// // Time goes from 0 to 5 seconds at 8000 Hz.
// stacksize('max');
// tfrsp(chirp([0:1/8000:5], 200, 2, 500, "logarithmic")',1:40001,128,'plot');
//  Authors
// 2001-08-31 Paul Kienzle pkienzle@users.sf.net

// * Fix documentation for quadratic case
// Copyright (C) 1999-2000 Paul Kienzle
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; If not, see <http://www.gnu.org/licenses/>.


  [nargout,nargin]=argn(0);
  if nargin < 1 | nargin > 6
    error("y = chirp(t [, f0 [, t1 [, f1 [, form [, phase]]]]])");
  end
  if nargin < 2, f0 = []; end
  if nargin < 3, t1 = []; end
  if nargin < 4, f1 = []; end
  if nargin < 5, form = []; end
  if nargin < 6, phase = []; end

  if isempty(f0), f0 = 0; end
  if isempty(t1), t1 = 1; end
  if isempty(f1), f1 = 100; end
  if isempty(form), form = "linear"; end
  if isempty(phase), phase = 0; end

  phase = 2*%pi*phase/360;

  if (form== "linear")
    a = %pi*(f1 - f0)/t1;
    b = 2*%pi*f0;
    y = cos(a*t.^2 + b*t + phase);
  elseif (form== "quadratic")
    a = (2/3*%pi*(f1-f0)/t1/t1);
    b = 2*%pi*f0;
    y = cos(a*t.^3 + b*t + phase);
  elseif (form== "logarithmic")
    a = 2*%pi*t1/log(f1-f0);
    b = 2*%pi*f0;
    x = (f1-f0)^(1/t1);
    y = cos(a*x.^t + b*t + phase);
  else
    error(sprintf("chirp doesnt understand ''%s''",form));
  end

endfunction

