/* Copyright (c) 2014 Benedikt Bitterli <benedikt.bitterli (at) gmail (dot) com>

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from
the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute
it freely, subject to the following restrictions:

    1. The origin of this software must not be misrepresented; you
       must not claim that you wrote the original software. If you
       use this software in a product, an acknowledgment in the
       product documentation would be appreciated but is not required.

    2. Altered source versions must be plainly marked as such, and
       must not be misrepresented as being the original software.

    3. This notice may not be removed or altered from any source
       distribution.

The Tungsten Renderer [online] [Accessed 2019]. 
Available from: https://github.com/tunabrain/tungsten
*/

/*
Copyright (c) 1998-2015, Matt Pharr, Greg Humphreys, and Wenzel Jakob.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this 
list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this 
list of conditions and the following disclaimer in the documentation and/or 
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

pbrt-v3 [online] [Accessed 2019]. 
Available from: https://github.com/mmp/pbrt-v3
*/

/* 
Copyrights (c) Side Effects Software Inc. All rights reserved.

Produced by:
 	Side Effects Software Inc
 	123 Front Street West, Suite 1401
 	Toronto, Ontario
 	Canada   M5J 2M2
 	416-504-9876

hair_eval.vfl [script, VEX] [Accessed 2019]
*/

#ifndef __extended_hair_h__
#define __extended_hair_h__

void assert(const int val) {}

// ----------------------------- ATTENUATION PART -----------------------------

float FrDielectric (float cosThetaI; float etaI; float etaT)
{
    float Rperp, Rparl;
    float cosTheta = clamp(cosThetaI, -1, 1);

    // the parallel polarized component amplitude
    float sinThetaI = 1.0 - cosTheta * cosTheta;

    float a = etaT * etaT - sinThetaI;
    if (a >= 0.0)
    {
        a = sqrt(a)/(etaT*cosTheta);
        Rparl = (a-etaI)/(a+etaI);
    } 
    else 
    {
        Rparl = 1.0;
    }

    // the perpendicular polarized one
    a = etaT * etaT - sinThetaI;
    if (a > 0.0)
    {
        a = etaT*cosTheta/sqrt(a);
        Rperp = (a-etaT)/(a+etaT);
    } 
    else
    {
        Rperp = 1.0;
    }
    float result = clamp((Rperp * Rperp + Rparl * Rparl) / 2, 0, 1);
    return result;
}

vector Ap(const float cosThetaO; const float eta; const float h; vector T; 
          const int p)
{
    float cosGammaO = sqrt(1 - h*h);
    float cosTheta = cosThetaO * cosGammaO;
    float f = FrDielectric(cosTheta, 1, eta);
    
    vector ap[];
    ap[0] = vector(f);
    ap[1] = (1-f)*(1-f) * T;
    ap[2] = ap[1] * T * f;

    return ap[p];
}

float ComputeApPdf(const float cosThetaO; const float eta; const float h; 
                   const vector sigma_a; const int p)
{
    float sinThetaO = sqrt(1 - cosThetaO*cosThetaO);
    float sinThetaT = sinThetaO/eta;
    float cosThetaT = sqrt(1 - sinThetaT*sinThetaT);
    
    float etap = sqrt(eta*eta - sinThetaO * sinThetaO) / cosThetaO;
    float sinGammaT = h / etap;
    float cosGammaT = sqrt(1 - sinGammaT * sinGammaT);
    
    vector T = exp(-sigma_a * (2 * cosGammaT / cosThetaT));
    
    vector ap[];
    float sumY = 0;
    for (int i = 0; i < 3; i++)
    {
        vector ap_temp = Ap(cosThetaO, eta, h, T, i);
        push(ap, ap_temp); 
        sumY += ap_temp.y;
    } 
    
    // normalizarion
    float apPdf = ap[p].y / sumY;

    return apPdf;
}

// ---------------------------- LONGITUDINAL PART -----------------------------

// Modified Bessel function of the first kind
float IO (const float x)
{
    float result = 1;
    float xSq = x*x;
    float xi = xSq;
    float denom = 4;
    for (int i = 1; i <= 10; ++i) 
    {
        result += (xi/denom);
        xi *= xSq;
        denom *= 4*(i + 1)*(i + 1);
    }
    return result;
}

float logIO (const float x)
{
    if (x > 12)
        return x + 0.5 * (log(1/(2*PI*x)) + 1/(8*x));
    else
        return log(IO(x));
}


float Mp(const float cosThetaI; const float cosThetaO; const float sinThetaI; 
         const float sinThetaO; const float l_variance)
{
    float a = cosThetaI * cosThetaO / l_variance;
    float b = sinThetaI * sinThetaO / l_variance;
    float mp = 0;
    if (l_variance <= 0.1)
    {   
        mp = exp(logIO(a) - b - 1 / l_variance + 0.6931 + log(1 / 
                (2 * l_variance))); 
    } 
    else
    {
        mp = (exp(-b) * IO(a)) / (sinh(1/l_variance) * 2 * l_variance);
    }
    return mp;
}

float SampleMp(const float cosThetaI; const float sinThetaI; const float sx; 
               const float sy; const float l_variance)
{
    float cosTheta = 1 + l_variance *log(sx + (1 - sx) * exp(-2/l_variance));
    float sinTheta = sqrt(1 - cosTheta*cosTheta);
    float cosPhi = cos(2 * PI * sy);
    return -cosTheta*sinThetaI + sinTheta*cosPhi*cosThetaI;
}

// ---------------------------- AZIMUTHAL PART --------------------------------

// "Light Scattering from Human Hair Fibers" 
float Phi(const int p; const float gammaO; const float gammaT)
{
    return 2 * p * gammaT - 2 * gammaO + p * PI;
}

// Logistic distribution is normalized and it is its own PDF
float Logistic(const float x; const float a_scale)
{
    float t_x = abs(x);
    return exp(-t_x / a_scale) / (a_scale * (1 + exp(-x / a_scale)) * 
           (1 + exp(-x / a_scale)));
}

float LogisticCDF(const float x; const float a_scale)
{
    return 1 / (1 + exp(-x/a_scale));
}

float TrimmedLogistic(const float x; const float a_scale; const float a; 
                      const float b)
{
    if (a < b)
    {
        return Logistic(x,a_scale) / (LogisticCDF(b,a_scale) - 
               LogisticCDF(a,a_scale));
    }
    return 0;
}

float sampleTrimmedLogistic(const float u; const float a_scale; const float a; 
                            const float b)
{
    if (a < b)
    {
        float k = LogisticCDF(b, a_scale) - LogisticCDF(a, a_scale);
        float x = -a_scale * log(1/(u * k + LogisticCDF(a, a_scale)) - 1);
        float x_t = clamp(x, a, b);
        return x_t;
    }
    return 0;
}

float Np(const float phi; const int p; const float a_scale; const float gammaO;
         const float gammaT)
{
    float dphi = phi - Phi(p, gammaO, gammaT);
    while (dphi > PI)
    {
        dphi -= 2*PI;
    }
    while (dphi < -PI)
    {
        dphi += 2 * PI;
    }
    return TrimmedLogistic(dphi, a_scale, -PI, PI);
}

// ------------------------------- ALBEDO PART --------------------------------

vector2 evaluateErfTerm(const float theta; const float mean; const float sdev)
{
    float	max_sdev = 10; 
    float	u, v;

    u = (2.0 * mean - theta) * (2.0 / (3.0 * PI));
    u *= 0.5;
    u += 0.5;

    v = sdev / max_sdev;

    return (vector2)vector(rawcolormap("erftable.rat", u, v,
				       "wrap", "clamp"));
}

// Albedo was computer to gaussian distributions which is close to logistic.
// ti is the incoming longitudinal angle.
float albedoGaussian (
    const float ti;
    const float lon_mean;
    const float lon_sdev;
    const float azi_mean;
    const float azi_sdev
    )
{
    if (abs(lon_mean) > PI / 2.0) return 0;
    if (lon_sdev < 0.0) return 0;
    if (abs(azi_mean) > PI) return 0;
    if (azi_sdev < 0.0) return 0;


    // f(ti, ai, tr, ar, u, v, p, q) :=
    //   exp(-((ti + tr) / 2 - u)^2 / (2 * v^2)) *
    //   exp(-(q - (ar - ai))^2 / (2 * p^2));
    //
    // g(ti, ai, u, v) :=
    //   integrate(
    //     integrate(
    //       f(ti, ai, tr, ar, u, v, p) * cos(tr),
    //       tr, -%pi / 2, %pi / 2),
    //     ar, ai, ai + 2 * %pi);
    //
 
    vector2	 value;
    float	 k1, k2, diff;

    k1 = erf(azi_mean / (sqrt(2.0) * azi_sdev)) -
	erf((azi_mean - 2.0 * PI) / (sqrt(2.0) * azi_sdev));

    k2 = 2.0 * lon_mean - ti;

    value = evaluateErfTerm(ti, lon_mean, lon_sdev);

    diff = 2.0 * value.x;

    diff *= PI * lon_sdev * azi_sdev;
    diff = max(0.0, diff);

    return diff / (2.0 * PI);
}

// Albedo for R lobe, which do not have inner intersections.
// Independent to the incoming azimuthal angle.
float
albedoCos(const float ti; const float mean; const float sdev)
{
    if (abs(mean) > PI / 2.0) return 0;
    if (sdev < 0.0) return 0;

    // f(ti, tr, ar, u, v) :=
    //   exp(-((ti + tr) / 2 - u)^2 / (2 * v^2)) * cos(ar / 2);
    //
    // g(ti, u, v) :=
    //   2 * integrate(
    //     integrate(
    //       f(ti, tr, ar, u, v) * cos(tr),
    //       tr, -%pi / 2, %pi / 2),
    //     ar, 0, %pi);
    //

    vector2	 value;
    float	 k, diff;

    k = 2.0 * mean - ti;
    value = evaluateErfTerm(ti, mean, sdev);

    diff = 2.0 * value.x;

    diff *= 2.0 * sqrt(2*PI) * sdev;
    diff = max(0.0, diff);

    return diff / (2.0 * PI);
}

#endif
