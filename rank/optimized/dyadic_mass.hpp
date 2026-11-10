// One-sided binary-grid bounds for the fourteen-fours mass inequalities.
struct DyadicMass {
    static constexpr unsigned scale = 1U << 20;
    unsigned outside_lower = 0;
    unsigned inside_upper = 0;
    unsigned sixes = 0;
    unsigned eights = 0;

    static unsigned reciprocal_upper(unsigned denominator) {
        if (denominator > 1024) return 1;
        const unsigned square = denominator * denominator;
        return (scale + square - 1) / square;
    }

    void append(unsigned b) {
        // The calling kernels clamp large denominators to 1025. Their
        // lower contribution is zero, and one upper unit covers either residue.
        sixes += b == 6;
        eights += b == 8;
        if (b > 1024) {
            ++inside_upper;
        } else if (b % 4 != 0) {
            outside_lower += scale / (b*b);
        } else if (b != 4) {
            inside_upper += reciprocal_upper(b);
        }
    }

    bool possible(unsigned rank, unsigned fours, unsigned last, unsigned slots,
                  unsigned root_lower, bool height_possible, bool only_235) const {
        if (fours >= 9 && !height_possible) return false;
        if (fours != 14) return true;
        if (rank <= 31 && (!only_235 || sixes > 1 || eights + slots < 5 || (last > 8 && eights < 5)))
            return false;
        if (last <= 4) return true;
        if (outside_lower >= scale/16) return false;
        // The next multiple of four cannot be smaller than this value.
        const unsigned first_multiple = 4*((last + 3)/4);
        const unsigned tail = slots * reciprocal_upper(first_multiple);
        // The nonunit prefix lcm divides the eventual root M; at a completed
        // pair root_lower is M itself (either bound may be clamped to 1025).
        const unsigned vacuum = reciprocal_upper(root_lower);
        return inside_upper + tail + vacuum > scale/16;
    }
};
