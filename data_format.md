## 2.1 The Data

The development of the customized fishplot function is based on the Morita et al. (2020) dataset.
This dataset includes 153 acute myeloid leukemia (AML) point mutation trees from 123 patients,
obtained using scDNA-seq. Tree inference was performed using the SCITE method.

The data are organized as a nested JSON-like structure representing clonal evolution trees with
associated metadata.

### Sample naming

Each sample is named using:
AML-{PATIENT_ID}-{TIMEPOINT}

Example: AML-63-005

### Tree structure

- Each sample corresponds to a single tree
- Nodes are hierarchically linked
- Each node may have zero or more children
- The root node represents the ancestral clone

### Node attributes

| Attribute        | Type              | Description |
|------------------|-------------------|-------------|
| node_id          | string / integer  | Unique node identifier |
| label            | list (optional)   | Typically ["Root"] |
| name             | string            | String representation |
| matching_label   | integer           | Mapping identifier |
| size_percent     | float (optional)  | Proportion of cells |
| gene_events      | dictionary        | Gene → mutation mapping |

### Example JSON

```json
{
  "AML-63-005": {
    "tree": {
      "node_id": "4",
      "label": ["Root"],
      "children": [
        {
          "node_id": "1",
          "size_percent": 0.511,
          "gene_events": {
            "IDH2": { "SNV": "p.R140Q" }
          }
        }
      ]
    },
    "metadata": {
      "VitalStatus": "Alive NOS"
    }
  }
}
```


