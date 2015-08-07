{
    "jobs": {
        "split_reference*": {
            "time": "00:10:00",
            "priority": "highprio"
        },
        "process_junctions":{
            "time": "03:00:00"
        },
        "augustus_*":{
            "time": "20h"
        },
        "merge*": {
            "time": "00:10:00"
        },
        "geneid_*":{
            "time": "3h", 
            "threads": "2"
        },
        "Genemark":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "48h",
            "priority": "lowprio"
        },
        "Genemark-ET":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "48h",
            "priority": "lowprio"
        },
        "Glimmer":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "48h",
            "priority": "lowprio"
        },
        "spaln":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "20h"
        },
        "Clean_transcripts": {
            "time": "00:10:00",
            "priority": "highprio"
        },
        "Pasa":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "72h",
            "priority": "lowprio"
        },
        "Transdecoder":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "12h"
        },
        "EVM_*":{
            "queue": "himem",
            "time": "90h",
            "priority": "lowprio",
            "threads": 16
        },
        "Select_evm":{
            "time": "3h"   
        },
        "Annotation_Update*":{
            "time": "24h",
            "threads": 2
        },
        "cmsearch":{
            "queue": "himem",
            "time": "48h",
            "priority": "lowprio"
        },
        "ncRNA_*":{
             "time": "00:10:00"
        },
        "Blast*":{
            "time": "3h"
        },
        "lncRNA":{
            "time": "1h"
        }
    } 
}

