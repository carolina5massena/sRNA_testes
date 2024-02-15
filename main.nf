process sayHello {
    input: 
        val x
    output:
        stdout
    
    script:
    """
    echo $x
    """
}

workflow {
  Channel.of('Bonjour', 'Ciao', 'Hello', 'Hola') | sayHello | view
}