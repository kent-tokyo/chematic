"""Independent formula token observer; D/T aliases retain isotope meaning."""
from collections import Counter
import re

TOKEN=re.compile(r'(?:\[(\d+)([A-Z][a-z]?|\*)\]|([A-Z][a-z]?|\*))(\d*)')


def observe(text):
    if not isinstance(text,str) or not text or text!=text.strip():
        raise ValueError('missing or noncanonical formula')
    charge=0;match=re.search(r'(?:(\+{2,}|-{2,})|([+-])(\d*))$',text)
    if match:
        if match[1]:charge=(1 if match[1][0]=='+' else -1)*len(match[1])
        else:charge=(1 if match[2]=='+' else -1)*int(match[3] or 1)
        text=text[:match.start()]
    if text=='H0':return [],charge
    counts=Counter();position=0
    while position<len(text):
        match=TOKEN.match(text,position)
        if not match:raise ValueError('unreadable formula token')
        mass=int(match[1]) if match[1] else None
        if mass is not None and not 0<mass<=65535:raise ValueError('unsupported isotope label')
        symbol=match[2] or match[3]
        if symbol in ('D','T'):
            if mass is not None:raise ValueError('isotope of isotope alias')
            symbol,mass='H',2 if symbol=='D' else 3
        count=int(match[4] or 1)
        if count<=0:raise ValueError('nonpositive count')
        counts[(symbol,mass)]+=count;position=match.end()
    if not counts:raise ValueError('empty formula')
    return [dict(element=s,isotope=m,count=n) for (s,m),n in sorted(counts.items(),key=lambda kv:(kv[0][0],kv[0][1] or 0))],charge


def totals(counts):
    result=Counter()
    for item in counts:result[item['element']]+=item['count']
    return result


def outcome(case,result):
    if type(result) is not dict or type(result.get('ok')) is not bool:raise ValueError('invalid envelope')
    if not result['ok']:
        if 'formula' in result or not isinstance(result.get('error'),str):raise ValueError('partial refusal')
        return dict(status='refused')
    try:counts,charge=observe(result['formula'])
    except (ValueError,KeyError,TypeError) as error:return dict(status='unreadable',error=str(error))
    same_counts=counts==case['expected_counts'];same_charge=charge==case['expected_charge']
    return dict(status='composition_match' if same_counts and same_charge else 'changed',
                isotope_counts_match=same_counts,charge_match=same_charge,
                element_totals_match=totals(counts)==totals(case['expected_counts']))
