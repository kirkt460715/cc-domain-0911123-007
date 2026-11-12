import json
from pathlib import Path

def load():
    return json.loads((Path(__file__).with_name("context.json")).read_text())

def run():
    from metrology.service import Service
    from metrology.ingress.readers import read_jsonl
    data = load()
    with Service() as app:
        app.catalog.add_source('bench')
        values = []
        for method in ['reference','portable']:
            app.catalog.add_instrument({'id':method,'label':'合成测量方法','channels':{'p':'kPa'}})
            app.catalog.add_certificate({'id':'cert-'+method,'instrument':method,'channel':'p','valid_from':'2026-01-01T00:00:00Z','valid_to':'2027-01-01T00:00:00Z','input_unit':'kPa','output_unit':'kPa','slope':'1','intercept':'0','uncertainty':'.01','minimum':'0','maximum':'100','issued_at':'2025-12-01T00:00:00Z'})
        app.catalog.add_policy({'id':'wide','unit':'kPa','lower':'0','upper':'100'})
        for method in ['reference','portable']:
            app.catalog.bind(method,'p','wide')
        for i,row in enumerate(data['records']):
            for j,method in enumerate(['reference','portable']):
                values.append({'id':f'{method}-{i}','source':'bench','sequence':i*2+j,'instrument':method,'channel':'p','event_time':row['at'],'value':row[method],'unit':data['unit'],'tags':{'specimen':row['specimen'],'repeat':str(row['repeat'])}})
        result = app.pipeline.ingest('paired', read_jsonl('\n'.join(json.dumps(r) for r in values)), received_at='2026-11-12T00:00:00Z')
        assert result['accepted'] == len(values)
        return {'ok': app.inspector.check()['invariants'] == 'ok', 'records': result['accepted']}


if __name__ == "__main__":
    print(json.dumps(run(), ensure_ascii=False, sort_keys=True))
