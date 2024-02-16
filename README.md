# sRNA_testes

# Dependências

- [Docker](https://www.docker.com/)
- [Nextflow.io](https://www.nextflow.io/)

# Pipeline

<div class="center">
   <img src="./docs/dag.png" style="text-align: center; width: 70%; border: 1px;margin: auto"/>
</dib>

# Run 

```bash
nextflow run -resume main.nf -with-report outputs/report.html -with-timeline outputs/timeline.html -with-dag docs/dag.png
```