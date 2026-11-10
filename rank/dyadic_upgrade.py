"""Reproducibly add the second-generation dyadic filters to both frozen kernels."""
from pathlib import Path


def one(text, old, new):
    if text.count(old) != 1:
        raise ValueError('incompatible frozen source: '+old)
    return text.replace(old,new)


def strengthened(name, source, header):
    text = one(source, 'struct SquarePrefix {', header+'\nstruct SquarePrefix {\n    DyadicMass dyadic;')
    text = one(text, 'SquarePrefix next=*this;', 'SquarePrefix next=*this; next.dyadic.append(b);')
    text = one(text, 'bool possible(unsigned rank,unsigned last,unsigned slots) const {',
               'bool possible(unsigned rank,unsigned last,unsigned slots,unsigned root_lower,bool height_possible,bool only_235) const {\n'
               '        if (!dyadic.possible(rank,fours,last,slots,root_lower,height_possible,only_235)) return false;')
    if name == 'fast':
        text = one(text, 'end.possible(rank_limit,last,0)',
                   'end.possible(rank_limit,last,0,static_cast<unsigned>(std::min<Wide>(root,1025)),root%8==0,(real_support & ~Mask(7))==0)')
        return one(text, 'square.possible(rank_limit,std::min<U>(prefix.back(),1025),h)',
                   'square.possible(rank_limit,std::min<U>(prefix.back(),1025),h,std::min<U>(L,1025),h>2 || L%8==0,(support & ~Mask(7))==0)')
    if name == 'exact':
        text = one(text, 'end.possible(rank,last,0)',
                   'end.possible(rank,last,0,s>1024?1025:s.get_ui(),divides(8,s),(sm & ~Mask(7))==0)')
        return one(text, 'state.possible(rank,prefix.back()>1024?1025:prefix.back().get_ui(),h)',
                   'state.possible(rank,prefix.back()>1024?1025:prefix.back().get_ui(),h,'
                   'L>1024?1025:L.get_ui(),h>2 || divides(8,L),(support & ~Mask(7))==0)')
    raise ValueError('unknown kernel')
