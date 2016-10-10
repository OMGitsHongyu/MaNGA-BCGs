test:
	for number in {0..79}; do \
		echo $$number; \
		python check_bcg_xyvviar_mpl5.py $$number; \
	done
