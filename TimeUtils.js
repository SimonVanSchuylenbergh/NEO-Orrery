export function MJDToJ2000(mjd) {
    return mjd - 51544.5;
}

export function J2000ToMJD(j2000) {
    return j2000 + 51544.5;
}

export function JDToMJD(JD) {
    return JD - 2400000.5;
}

export function MJDToJD(MJD) {
    return MJD + 2400000.5;
}

export function datetimeToMJD(datetimeStr) {
    // Split the date and time parts
    const [datePart, timePart] = datetimeStr.split(' ');

    // Extract year, month, and day from the date part
    const [year, month, day] = datePart.split('-').map(Number);

    // Extract hours, minutes, and seconds from the time part
    const [hours, minutes, seconds] = timePart.split(':').map(Number);

    // Adjust month and year for January and February
    let y = year;
    let m = month;
    if (m <= 2) {
        y -= 1;
        m += 12;
    }

    // Calculate Julian Date (JD)
    const A = Math.floor(y / 100);
    const B = 2 - A + Math.floor(A / 4);
    const JD = Math.floor(365.25 * (y + 4716)) + Math.floor(30.6001 * (m + 1)) + day + B - 1524.5
        + (hours + minutes / 60 + seconds / 3600) / 24;

    // Convert to Modified Julian Date (MJD)
    const MJD = JD - 2400000.5;

    return MJD;
}


export function MJDToDatetime(mjd) {
    // Convert MJD to Julian Date (JD)
    const JD = mjd + 2400000.5;

    // Convert JD to Gregorian date
    const F = JD + 0.5;
    const Z = Math.floor(F);
    const frac = F - Z;
    
    let A = Z;
    if (Z >= 2299161) {
        const alpha = Math.floor((Z - 1867216.25) / 36524.25);
        A = Z + 1 + alpha - Math.floor(alpha / 4);
    }

    const B = A + 1524;
    const C = Math.floor((B - 122.1) / 365.25);
    const D = Math.floor(365.25 * C);
    const E = Math.floor((B - D) / 30.6001);

    const day = B - D - Math.floor(30.6001 * E);
    let month = (E < 14) ? E - 1 : E - 13;
    let year = (month > 2) ? C - 4716 : C - 4715;

    // Calculate hours, minutes, and seconds from the fractional day
    const dayFrac = frac * 24;
    const hours = Math.floor(dayFrac);
    const minuteFrac = (dayFrac - hours) * 60;
    const minutes = Math.floor(minuteFrac);
    const seconds = Math.floor((minuteFrac - minutes) * 60);

    // Zero-padding for hours, minutes, and seconds
    const pad = (num) => num.toString().padStart(2, '0');

    return `${year}-${pad(month)}-${pad(day)} ${pad(hours)}:${pad(minutes)}:${pad(seconds)}`;
}


export function datetimeToJ2000(datetime) {
    return MJDToJ2000(datetimeToMJD(datetime));
}

export function J2000ToDatetime(j2000) {
    return MJDToDatetime(J2000ToMJD(j2000))
}


