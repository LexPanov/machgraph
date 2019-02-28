float sdSphere(vec3 p, float s) {
    return length(p) - s;
}

vec2 opU(vec2 d1, vec2 d2) {
    return (d1.x < d2.x) ? d1 : d2;
}

#define ZERO (min(iFrame,0))

vec2 map(in vec3 pos) {
    vec2 res = vec2(1e10, 0.0);

    res = opU(res, vec2(sdSphere(pos - vec3(-0.6, 0.25, 0.75), 0.25), 46.9));

    return res;
}

// http://iquilezles.org/www/articles/boxfunctions/boxfunctions.htm
vec2 iBox(in vec3 ro, in vec3 rd, in vec3 rad) {
    vec3
    m = 1.0 / rd;
    vec3
    n = m * ro;
    vec3
    k = abs(m) * rad;
    vec3
    t1 = -n - k;
    vec3
    t2 = -n + k;
    return vec2(max(max(t1.x, t1.y), t1.z),
                min(min(t2.x, t2.y), t2.z));
}

const float maxHei = 0.8;

vec2 castRay(in vec3 ro, in vec3 rd) {
    vec2 res = vec2(-1.0, -1.0);

    float tmin = 1.0;
    float tmax = 20.0;

    // raytrace floor plane
    float tp1 = (0.0 - ro.y) / rd.y;
    if (tp1 > 0.0) {
        tmax = min(tmax, tp1);
        res = vec2(tp1, 1.0);
    }
    //else return res;

    // raymarch primitives
    vec2 tb = iBox(ro - vec3(0.0, 0.4, 0.0), rd, vec3(2.5, 0.41, 2.5));
    if (tb.x < tb.y && tb.y > 0.0 && tb.x < tmax) {
        tmin = max(tb.x, tmin);
        tmax = min(tb.y, tmax);

        float t = tmin;
        for (int i = 0; i < 70 && t < tmax; i++) {
            vec2 h = map(ro + rd * t);
            if (abs(h.x) < (0.0001 * t)) {
                res = vec2(t, h.y);
                break;
            }
            t += h.x;
        }
    }

    return res;
}


// http://iquilezles.org/www/articles/rmshadows/rmshadows.htm
float calcSoftshadow(in vec3 ro, in vec3 rd, in float mint, in float tmax) {
    // bounding volume
    float tp = (maxHei - ro.y) / rd.y;
    if (tp > 0.0) tmax = min(tmax, tp);

    float res = 1.0;
    float t = mint;
    for (int i = ZERO; i < 16; i++) {
        float h = map(ro + rd * t).x;
        res = min(res, 8.0 * h / t);
        t += clamp(h, 0.02, 0.10);
        if (res < 0.005 || t > tmax) break;
    }
    return clamp(res, 0.0, 1.0);
}

// http://iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
vec3 calcNormal(in vec3 pos) {

    vec2 e = vec2(1.0, -1.0) * 0.5773 * 0.0005;
    return normalize(e.xyy * map(pos + e.xyy).x +
                     e.yyx * map(pos + e.yyx).x +
                     e.yxy * map(pos + e.yxy).x +
                     e.xxx * map(pos + e.xxx).x);

}

float calcAO(in vec3 pos, in vec3 nor) {
    float occ = 0.0;
    float sca = 1.0;
    for (int i = ZERO; i < 5; i++) {
        float hr = 0.01 + 0.12 * float(i) / 4.0;
        vec3
        aopos = nor * hr + pos;
        float dd = map(aopos).x;
        occ += -(dd - hr) * sca;
        sca *= 0.95;
    }
    return clamp(1.0 - 3.0 * occ, 0.0, 1.0) * (0.5 + 0.5 * nor.y);
}

// http://iquilezles.org/www/articles/checkerfiltering/checkerfiltering.htm
float checkersGradBox(in vec2 p) {
    // filter kernel
    vec2
    w = fwidth(p) + 0.001;
    // analytical integral (box filter)
    vec2
    i = 2.0 * (abs(fract((p - 0.5 * w) * 0.5) - 0.5) - abs(fract((p + 0.5 * w) * 0.5) - 0.5)) / w;
    // xor pattern
    return 0.5 - 0.5 * i.x * i.y;
}

vec3 render(in vec3 ro, in vec3 rd) {
    vec3 col = vec3(0.7, 0.9, 1.0) + rd.y * 0.8;
    vec2 res = castRay(ro, rd);
    float t = res.x;
    float m = res.y;
    if (m > -0.5) {
        vec3
        pos = ro + t * rd;
        vec3
        nor = (m < 1.5) ? vec3(0.0, 1.0, 0.0) : calcNormal(pos);
        vec3
        ref = reflect(rd, nor);

        // material
        col = 0.45 + 0.35 * sin(vec3(0.05, 0.08, 0.10) * (m - 1.0));
        if (m < 1.5) {

            float f = checkersGradBox(5.0 * pos.xz);
            col = 0.3 + f * vec3(0.1);
        }

        // lighting
        float occ = calcAO(pos, nor);
        vec3 lig = normalize(vec3(-0.4, 0.7, -0.6));
        vec3 hal = normalize(lig - rd);
        float amb = clamp(0.5 + 0.5 * nor.y, 0.0, 1.0);
        float dif = clamp(dot(nor, lig), 0.0, 1.0);
        float bac = clamp(dot(nor, normalize(vec3(-lig.x, 0.0, -lig.z))), 0.0, 1.0) * clamp(1.0 - pos.y, 0.0, 1.0);
        float dom = smoothstep(-0.2, 0.2, ref.y);
        float fre = pow(clamp(1.0 + dot(nor, rd), 0.0, 1.0), 2.0);

        dif *= calcSoftshadow(pos, lig, 0.02, 2.5);
        dom *= calcSoftshadow(pos, ref, 0.02, 2.5);

        float spe = pow(clamp(dot(nor, hal), 0.0, 1.0), 16.0) *
                    dif *
                    (0.04 + 0.96 * pow(clamp(1.0 + dot(hal, rd), 0.0, 1.0), 5.0));

        vec3 lin = vec3(0.0);
        lin += 1.40 * dif * vec3(1.00, 0.80, 0.55);
        lin += 0.20 * amb * vec3(0.40, 0.60, 1.00) * occ;
        lin += 0.40 * dom * vec3(0.40, 0.60, 1.00) * occ;
        lin += 0.50 * bac * vec3(0.25, 0.25, 0.25) * occ;
        lin += 0.25 * fre * vec3(1.00, 1.00, 1.00) * occ;
        col = col * lin;
        col += 9.00 * spe * vec3(1.00, 0.90, 0.70);

        col = mix(col, vec3(0.8, 0.9, 1.0), 1.0 - exp(-0.0002 * t * t * t));
    }

    return vec3(clamp(col, 0.0, 1.0));
}

mat3 setCamera(in vec3 ro, in vec3 ta, float cr) {
    vec3 cw = normalize(ta - ro);
    vec3 cp = vec3(sin(cr), cos(cr), 0.0);
    vec3 cu = normalize(cross(cw, cp));
    vec3 cv = (cross(cu, cw));
    return mat3(cu, cv, cw);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2
    mo = iMouse.xy / iResolution.xy;
    float time = 15.0 + iTime;

    // camera
    vec3 ro = vec3(4.6 * cos(0.1 * time + 6.0 * mo.x), 1.0 + 2.0 * mo.y, 0.5 + 4.6 * sin(0.1 * time + 6.0 * mo.x));
    vec3 ta = vec3(-0.5, -0.2, 0.8);
    // camera-to-world transformation
    mat3 ca = setCamera(ro, ta, 0.0);

    vec3 tot = vec3(0.0);

    vec2
    p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;

    // ray direction
    vec3 rd = ca * normalize(vec3(p.xy, 4.0));

    // render
    vec3 col = render(ro, rd);

    // gamma
    col = pow(col, vec3(0.4545));

    tot += col;


    fragColor = vec4(tot, 1.0);
}