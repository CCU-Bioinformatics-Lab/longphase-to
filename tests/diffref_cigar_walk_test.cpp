// Checks getWindowsDiffRef (ParsingBam.cpp): for a variant on a read, it records the offset of every base
// within 100 aligned bases on each side that differs from the reference. Offset k means the k-th aligned
// base (M, = or X) away from the variant; insertions, deletions and clips take no step.
// The same alignment must give the same offsets whether its CIGAR is written with M or with =/X.
// Built and run by tests/test_diffref_cigar_walk.sh, linked against the project's object files.
#include <htslib/sam.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

void getWindowsDiffRef(const uint32_t *cigar, int cigarIndex, const bam1_t &aln, const std::string &refString,
                       int readPos, int readOffset, int refPos, std::vector<int8_t> &offsets);

namespace {

uint32_t lcgState = 12345u;
int nextInt(int n) {
    lcgState = lcgState * 1103515245u + 12345u;
    return static_cast<int>((lcgState >> 8) % static_cast<uint32_t>(n));
}

struct Read {
    std::string seq;
    std::vector<std::pair<int, int> > ops;  // (BAM op, length) with M for aligned bases
};

void addOp(Read &read, int op) {
    if (!read.ops.empty() && read.ops.back().first == op) {
        read.ops.back().second++;
    } else {
        read.ops.push_back(std::make_pair(op, 1));
    }
}

// Read over ref: mismatches, and deletions / insertions (each with the given chance in percent).
Read makeRead(const std::string &ref, int delPercent, int insPercent) {
    const char *bases = "ACGT";
    Read read;
    int clip = nextInt(3);
    for (int k = 0; k < clip; k++) { read.seq += bases[nextInt(4)]; addOp(read, BAM_CSOFT_CLIP); }
    for (int r = 0; r < static_cast<int>(ref.size());) {
        int u = nextInt(100);
        if (u < delPercent && r > 5 && r < static_cast<int>(ref.size()) - 5) {
            int len = 1 + nextInt(4);
            for (int k = 0; k < len; k++) addOp(read, BAM_CDEL);
            r += len;
            continue;
        }
        if (u < delPercent + insPercent && r > 5 && r < static_cast<int>(ref.size()) - 5) {
            int len = 1 + nextInt(4);
            for (int k = 0; k < len; k++) { read.seq += bases[nextInt(4)]; addOp(read, BAM_CINS); }
        }
        char b = ref[r];
        if (nextInt(100) < 8) b = bases[(std::string(bases).find(b) + 1 + nextInt(3)) % 4];
        read.seq += b;
        addOp(read, BAM_CMATCH);
        r++;
    }
    clip = nextInt(3);
    for (int k = 0; k < clip; k++) { read.seq += bases[nextInt(4)]; addOp(read, BAM_CSOFT_CLIP); }
    return read;
}

void pushCigar(std::vector<uint32_t> &cigar, int op, uint32_t len) {
    if (!cigar.empty() && static_cast<int>(bam_cigar_op(cigar.back())) == op) {
        cigar.back() += len << BAM_CIGAR_SHIFT;
    } else {
        cigar.push_back(bam_cigar_gen(len, op));
    }
}

std::vector<uint32_t> cigarM(const Read &read) {
    std::vector<uint32_t> cigar;
    for (size_t k = 0; k < read.ops.size(); k++) pushCigar(cigar, read.ops[k].first, read.ops[k].second);
    return cigar;
}

std::vector<uint32_t> cigarEqx(const Read &read, const std::string &ref) {
    std::vector<uint32_t> cigar;
    int q = 0, t = 0;
    for (size_t k = 0; k < read.ops.size(); k++) {
        int op = read.ops[k].first, len = read.ops[k].second;
        if (op == BAM_CMATCH) {
            for (int j = 0; j < len; j++) pushCigar(cigar, read.seq[q + j] == ref[t + j] ? BAM_CEQUAL : BAM_CDIFF, 1);
            q += len; t += len;
        } else {
            pushCigar(cigar, op, len);
            if (op == BAM_CINS || op == BAM_CSOFT_CLIP) q += len;
            if (op == BAM_CDEL) t += len;
        }
    }
    return cigar;
}

// Calls getWindowsDiffRef the way the BAM parser does: on the operation that holds the variant.
std::vector<int> walk(const std::vector<uint32_t> &cigar, const std::string &ref, const std::string &seq, int variantPos) {
    bam1_t *aln = bam_init1();
    bam_set1(aln, 1, "r", 0, 0, 0, 60, cigar.size(), cigar.data(), -1, -1, 0, seq.size(), seq.c_str(), NULL, 0);
    int q = 0, t = 0, opIndex = -1, opReadStart = 0, offset = 0;
    for (size_t k = 0; k < cigar.size(); k++) {
        int op = bam_cigar_op(cigar[k]), len = bam_cigar_oplen(cigar[k]);
        bool aligned = op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF;
        if (aligned && variantPos >= t && variantPos < t + len) { opIndex = k; opReadStart = q; offset = variantPos - t; break; }
        if (aligned) { q += len; t += len; }
        else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) q += len;
        else if (op == BAM_CDEL) t += len;
    }
    std::vector<int8_t> offsets;
    getWindowsDiffRef(cigar.data(), opIndex, *aln, ref, opReadStart, offset, variantPos, offsets);
    bam_destroy1(aln);
    std::vector<int> result(offsets.begin(), offsets.end());
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<int> expected(const Read &read, const std::string &ref, int variantPos) {
    std::vector<std::pair<int, int> > aligned;  // (read position, reference position)
    int q = 0, t = 0;
    for (size_t k = 0; k < read.ops.size(); k++) {
        int op = read.ops[k].first, len = read.ops[k].second;
        if (op == BAM_CMATCH) { for (int j = 0; j < len; j++) aligned.push_back(std::make_pair(q + j, t + j)); q += len; t += len; }
        else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) q += len;
        else if (op == BAM_CDEL) t += len;
    }
    int at = -1;
    for (size_t k = 0; k < aligned.size(); k++) if (aligned[k].second == variantPos) { at = k; break; }
    std::vector<int> result;
    for (int d = 1; d <= 100; d++) {
        if (at - d >= 0 && read.seq[aligned[at - d].first] != ref[aligned[at - d].second]) result.push_back(-d);
        if (at + d < static_cast<int>(aligned.size()) && read.seq[aligned[at + d].first] != ref[aligned[at + d].second]) result.push_back(d);
    }
    std::sort(result.begin(), result.end());
    return result;
}

std::string show(const std::vector<int> &v) {
    std::string s;
    for (size_t k = 0; k < v.size(); k++) { char buf[16]; snprintf(buf, sizeof(buf), " %d", v[k]); s += buf; }
    return s;
}

// Returns the number of reads whose offsets differ from the expected ones, for M and for =/X CIGARs.
int runCase(const char *name, int delPercent, int insPercent, int reads) {
    int badM = 0, badEqx = 0;
    for (int n = 0; n < reads; n++) {
        std::string ref;
        for (int i = 0; i < 400; i++) ref += "ACGT"[nextInt(4)];
        Read read = makeRead(ref, delPercent, insPercent);
        std::vector<int> refPositions;
        { int t = 0; for (size_t k = 0; k < read.ops.size(); k++) {
              if (read.ops[k].first == BAM_CMATCH) for (int j = 0; j < read.ops[k].second; j++) refPositions.push_back(t + j);
              if (read.ops[k].first == BAM_CMATCH || read.ops[k].first == BAM_CDEL) t += read.ops[k].second; } }
        int variantPos = refPositions[120 + nextInt(static_cast<int>(refPositions.size()) - 240)];
        std::vector<int> want = expected(read, ref, variantPos);
        std::vector<int> gotM = walk(cigarM(read), ref, read.seq, variantPos);
        std::vector<int> gotEqx = walk(cigarEqx(read, ref), ref, read.seq, variantPos);
        if (gotM != want && ++badM == 1) printf("    %s, M example: expected%s\n      got%s\n", name, show(want).c_str(), show(gotM).c_str());
        if (gotEqx != want && ++badEqx == 1) printf("    %s, =/X example: expected%s\n      got%s\n", name, show(want).c_str(), show(gotEqx).c_str());
    }
    printf("%s: %s (%d reads; wrong offsets with M CIGAR: %d, with =/X CIGAR: %d)\n",
           (badM || badEqx) ? "FAIL" : "PASS", name, reads, badM, badEqx);
    return badM + badEqx;
}

}  // namespace

int main() {
    int bad = 0;
    bad += runCase("mismatches only", 0, 0, 300);
    bad += runCase("mismatches and deletions", 3, 0, 300);
    bad += runCase("mismatches and insertions", 0, 3, 300);
    bad += runCase("mismatches, deletions and insertions", 3, 3, 300);
    return bad ? 1 : 0;
}
