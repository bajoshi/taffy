for f in *.eps; do
  convert ./"$f" ./"${f%.eps}.png"
done

# to run you should cd into the folder where the script is and 
# . convert_eps_to_png.sh